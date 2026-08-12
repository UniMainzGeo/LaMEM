# check_valgrind_report.jl
#
# Scans all Valgrind XML reports (produced under `make grind`, one file per MPI rank, named
# <testname>_<pid>.xml) below the current directory, and reports:
#   1. "Leak_DefinitelyLost" errors (memory leaks)
#   2. "UninitValue" / "UninitCondition" errors (use of uninitialised values), including the
#      "origin" stack (where the value was created) that --track-origins=yes provides.
#
# Source file:line is shown for each stack frame where Valgrind has debug info. Occurrences are
# merged/deduplicated by the single source line that actually identifies the site - where the
# unfreed allocation was made (leaks), or where the uninitialised value is used (uninit errors) -
# rather than the full call stack, since that's what "the same leak/uninit site" means in practice,
# regardless of which (possibly varying) call path got there. This matters a lot, since the same site
# typically fires once per MPI rank and/or once per loop iteration, and would otherwise flood the
# report with near-identical copies.
#
# Exits with a non-zero status if anything is found, so `make grind` fails loudly and CI can catch it.
#
# Usage: julia --project=../. check_valgrind_report.jl

function find_xml_files(root::String)
    files = String[]
    for (dir, _, fnames) in walkdir(root)
        for f in fnames
            endswith(f, ".xml") && push!(files, joinpath(dir, f))
        end
    end
    return sort(files)
end

get_tag(block::AbstractString, tag::AbstractString) = begin
    m = match(Regex("<$(tag)>(.*?)</$(tag)>", "s"), block)
    isnothing(m) ? nothing : m.captures[1]
end

"""
    parse_frame(fblock::AbstractString)

Parse a single `<frame>...</frame>` block into a NamedTuple with function name, object (library/exe),
source file:line (when Valgrind has debug info for that frame), and the raw filename (Valgrind's
`<file>` tag is just the basename, e.g. "utils.c" - no directory), which lets us check it against
LaMEM's own source file list.
"""
function parse_frame(fblock::AbstractString)
    fn   = get_tag(fblock, "fn")
    obj  = get_tag(fblock, "obj")
    file = get_tag(fblock, "file")
    line = get_tag(fblock, "line")

    fn = isnothing(fn) ? "???" : fn

    loc = if !isnothing(file) && !isnothing(line)
        "$(file):$(line)"
    elseif !isnothing(obj)
        "in $(obj)"
    else
        "unknown location"
    end

    has_src = !isnothing(file) && !isnothing(line)

    return (fn=fn, loc=loc, has_src=has_src, file=file)
end

"""
    parse_stacks(block::AbstractString)

Extract all `<stack>...</stack>` sections in order (an error's own stack, then - for
--track-origins=yes uninitialised-value errors - the "origin" stack), each as a Vector of frames.
"""
function parse_stacks(block::AbstractString)
    stacks = Vector{Vector{NamedTuple}}()
    for sm in eachmatch(r"<stack>(.*?)</stack>"s, block)
        frames = [parse_frame(fm.captures[1]) for fm in eachmatch(r"<frame>(.*?)</frame>"s, sm.captures[1])]
        push!(stacks, frames)
    end
    return stacks
end

"""
    lamem_source_files(src_dir::String)

Collect the basenames of all LaMEM source/header files under `src_dir` (recursively). Stack frames are
then matched against this list of LaMEM's *actual* files, rather than guessing at library name
patterns (PETSc, MPI, BLAS, ...) - simpler and doesn't need updating if the dependency set changes.
"""
function lamem_source_files(src_dir::String)
    files = Set{String}()
    if !isdir(src_dir)
        return files
    end
    for (dir, _, fnames) in walkdir(src_dir)
        for f in fnames
            if any(endswith(f, ext) for ext in (".c", ".cpp", ".cc", ".h", ".hpp"))
                push!(files, f)
            end
        end
    end
    return files
end

is_lamem_frame(fr, lamem_files::Set{String}) = fr.has_src && !isnothing(fr.file) && (fr.file in lamem_files)

"""
    display_frames(frames::Vector, lamem_files::Set{String}; max_frames=6)

Pick which frames are worth showing/keying on. Frames inside LaMEM's own source (per `lamem_files`)
are prioritized over everything else - we're not trying to debug PETSc, MPI, BLAS/LAPACK, etc., even
when those ship with debug info and Valgrind happily reports a file:line inside them. Falls back to
any frame with known source info if no LaMEM frame is present, and finally to the raw frames if none
carry source info at all, so there's always something to look at.
"""
function display_frames(frames::Vector, lamem_files::Set{String}; max_frames=6)
    to_show = filter(f -> is_lamem_frame(f, lamem_files), frames)
    isempty(to_show) && (to_show = filter(f -> f.has_src, frames))  # fall back: any frame with source
    isempty(to_show) && (to_show = frames)                          # fall back: raw frames
    return first(to_show, max_frames)
end

function print_stack(frames::Vector; indent="      ")
    for (i, fr) in enumerate(frames)
        println("$(indent)#$(i-1) $(fr.fn) ($(fr.loc))")
    end
end

# A hashable signature for a (trimmed) stack, used to merge occurrences that are really "the same"
# leak/uninit site - e.g. hit once per MPI rank, or once per loop iteration.
frame_sig(frames::Vector) = Tuple("$(fr.fn) @ $(fr.loc)" for fr in frames)

"""
    primary_location(frames::Vector, lamem_files::Set{String})

The single source line that actually identifies a leak/uninit site: the first frame inside LaMEM's
own source (falling back to any frame with source info, then to the very first frame, if none of
LaMEM's own files appear in the stack at all). This - not the whole call stack - is what occurrences
get grouped/merged by, since the interesting thing for deduplication purposes is "where in LaMEM was
this unfreed allocation made" / "where in LaMEM is this uninitialised value used", regardless of which
(possibly varying) call path got there, and regardless of which library frames happen to sit above it.
"""
function primary_location(frames::Vector, lamem_files::Set{String})
    isempty(frames) && return (fn="???", loc="unknown location", has_src=false, file=nothing)
    idx = findfirst(f -> is_lamem_frame(f, lamem_files), frames)
    idx = isnothing(idx) ? findfirst(f -> f.has_src, frames) : idx
    return isnothing(idx) ? frames[1] : frames[idx]
end

function parse_leaks(content::AbstractString)
    leaks = NamedTuple[]
    for m in eachmatch(r"<error>(.*?)</error>"s, content)
        block = m.captures[1]
        get_tag(block, "kind") == "Leak_DefinitelyLost" || continue

        bytes_m  = match(r"<leakedbytes>(\d+)</leakedbytes>", block)
        blocks_m = match(r"<leakedblocks>(\d+)</leakedblocks>", block)
        bytes  = isnothing(bytes_m)  ? 0 : parse(Int, bytes_m.captures[1])
        blocks = isnothing(blocks_m) ? 0 : parse(Int, blocks_m.captures[1])

        stacks = parse_stacks(block)
        push!(leaks, (bytes=bytes, blocks=blocks, stack=isempty(stacks) ? [] : stacks[1]))
    end
    return leaks
end

function parse_uninit(content::AbstractString)
    issues = NamedTuple[]
    for m in eachmatch(r"<error>(.*?)</error>"s, content)
        block = m.captures[1]
        kind = get_tag(block, "kind")
        kind in ("UninitValue", "UninitCondition") || continue

        what = something(get_tag(block, "what"), "Use of uninitialised value")
        aux  = get_tag(block, "auxwhat")
        if isnothing(aux)
            aux = get_tag(block, "text")  # xauxwhat wraps its text in <text>...</text>
        end

        stacks = parse_stacks(block)
        error_stack  = isempty(stacks) ? [] : stacks[1]
        origin_stack = length(stacks) >= 2 ? stacks[2] : []

        push!(issues, (what=what, aux=aux, error_stack=error_stack, origin_stack=origin_stack))
    end
    return issues
end

"""
    format_files(files::Vector{String}; max_show=3)

Compact "seen in" listing: show a handful of example files, then "and N more" instead of dumping a
potentially huge list (one file per MPI rank per test that hit the same site).
"""
function format_files(files::Vector{String}; max_show=3)
    n = length(files)
    n <= max_show && return join(files, ", ")
    return join(files[1:max_show], ", ") * ", and $(n - max_show) more"
end

function group_leaks(entries, lamem_files::Set{String})
    # entries: Vector{Tuple{String,NamedTuple}} of (relative xml file, leak)
    groups = Dict{Any, Dict{Symbol,Any}}()
    for (file, leak) in entries
        frames = display_frames(leak.stack, lamem_files)
        key = frame_sig([primary_location(frames, lamem_files)])
        g = get!(groups, key) do
            Dict{Symbol,Any}(:count=>0, :bytes=>0, :blocks=>0, :frames=>frames, :files=>String[])
        end
        g[:count]  += 1
        g[:bytes]  += leak.bytes
        g[:blocks] += leak.blocks
        push!(g[:files], file)
    end
    return groups
end

function group_uninit(entries, lamem_files::Set{String})
    groups = Dict{Any, Dict{Symbol,Any}}()
    for (file, u) in entries
        err_frames    = display_frames(u.error_stack, lamem_files)
        origin_frames = display_frames(u.origin_stack, lamem_files)
        key = (u.what, frame_sig([primary_location(err_frames, lamem_files)]))
        g = get!(groups, key) do
            Dict{Symbol,Any}(:count=>0, :what=>u.what, :aux=>u.aux,
                              :error_frames=>err_frames, :origin_frames=>origin_frames, :files=>String[])
        end
        g[:count] += 1
        push!(g[:files], file)
    end
    return groups
end

function main()
    root = pwd()
    xmlfiles = find_xml_files(root)

    if isempty(xmlfiles)
        println("No Valgrind XML reports found under $(root).")
        println("Did you run 'make grind' first?")
        exit(1)
    end

    println("="^70)
    println("Valgrind report (scanned $(length(xmlfiles)) XML file(s) under $(root))")
    println("="^70)

    lamem_src_dir = joinpath(root, "..", "src")
    lamem_files   = lamem_source_files(lamem_src_dir)
    if isempty(lamem_files)
        println("\nWarning: no LaMEM source files (.c/.cpp/.cc/.h/.hpp) found under $(lamem_src_dir).")
        println("Stack frames won't be filtered down to LaMEM's own code - falling back to any frame with source info.")
    end

    all_leaks  = Tuple{String,NamedTuple}[]
    all_uninit = Tuple{String,NamedTuple}[]

    for file in xmlfiles
        content = read(file, String)
        rel = relpath(file, root)
        for leak in parse_leaks(content)
            push!(all_leaks, (rel, leak))
        end
        for u in parse_uninit(content)
            push!(all_uninit, (rel, u))
        end
    end

    leak_groups   = group_leaks(all_leaks, lamem_files)
    uninit_groups = group_uninit(all_uninit, lamem_files)

    if !isempty(leak_groups)
        println("\nMemory leaks (definitely lost) - $(length(leak_groups)) unique site(s), " *
                "$(length(all_leaks)) total occurrence(s):")
        for (_, g) in sort(collect(leak_groups), by = kv -> -kv[2][:count])
            println("\n  x$(g[:count])   $(g[:bytes]) bytes / $(g[:blocks]) blocks (total across occurrences)")
            print_stack(g[:frames])
            println("      seen in: $(format_files(unique(g[:files])))")
        end
    end

    if !isempty(uninit_groups)
        println("\nUninitialised-value errors - $(length(uninit_groups)) unique site(s), " *
                "$(length(all_uninit)) total occurrence(s):")
        for (_, g) in sort(collect(uninit_groups), by = kv -> -kv[2][:count])
            println("\n  x$(g[:count])   $(g[:what])")
            print_stack(g[:error_frames])
            if !isempty(g[:origin_frames])
                println("      origin: $(something(g[:aux], "created at"))")
                print_stack(g[:origin_frames]; indent="        ")
            end
            println("      seen in: $(format_files(unique(g[:files])))")
        end
    end

    println()
    println("="^70)
    any_failed = !isempty(all_leaks) || !isempty(all_uninit)
    if any_failed
        total_bytes = sum(leak.bytes for (_, leak) in all_leaks; init=0)
        println("FAILED:")
        !isempty(all_leaks) && println("  - $(length(all_leaks)) leak occurrence(s) " *
                "($(length(leak_groups)) unique site(s)), totalling $(total_bytes) bytes.")
        !isempty(all_uninit) && println("  - $(length(all_uninit)) uninitialised-value occurrence(s) " *
                "($(length(uninit_groups)) unique site(s)).")
        println("="^70)
        exit(1)
    else
        println("PASSED: no 'definitely lost' leaks or uninitialised-value errors found.")
        println("="^70)
        exit(0)
    end
end

main()
