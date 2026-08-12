#!/usr/bin/env julia
"""
petsc_getrestore_check.jl

Static checker that pairs PETSc *GetArray* / *RestoreArray* calls
(VecGetArray, VecGetArrayRead, DMDAVecGetArray, DMDAVecGetArrayRead, and
their DOF/Write variants) and flags any that are unmatched within a
function body.

Usage:
    julia petsc_getrestore_check.jl path/to/source/dir [more/dirs ...]

Notes / limitations:
  - This is a lightweight lexical/brace-scope scanner, not a full C parser.
    It ignores calls inside #if 0 / comments (comments are stripped first)
    but does not understand preprocessor branching (#ifdef), so a Get/Restore
    split across mutually-exclusive #ifdef branches will be reported as
    unmatched. Review flagged cases manually.
  - Matching key:
      VecGetArray*      -> the Vec argument (arg 0)
      DMDAVecGetArray*  -> the (DM, Vec) argument pair (args 0,1)
  - Matches are tracked as a stack per key per function, in source order,
    so nested/re-entrant Get calls on the same object are handled correctly.
  - Input is treated as ASCII source code: any non-ASCII byte is dropped
    before scanning (mirrors Python's `errors="ignore"` decoding), so all
    offsets stay simple 1-byte-per-character positions.
"""

const DOCSTRING = raw"""
petsc_getrestore_check.jl

Static checker that pairs PETSc *GetArray* / *RestoreArray* calls
(VecGetArray, VecGetArrayRead, DMDAVecGetArray, DMDAVecGetArrayRead, and
their DOF/Write variants) and flags any that are unmatched within a
function body.

Usage:
    julia petsc_getrestore_check.jl path/to/source/dir [more/dirs ...]

Notes / limitations:
  - This is a lightweight lexical/brace-scope scanner, not a full C parser.
    It ignores calls inside #if 0 / comments (comments are stripped first)
    but does not understand preprocessor branching (#ifdef), so a Get/Restore
    split across mutually-exclusive #ifdef branches will be reported as
    unmatched. Review flagged cases manually.
  - Matching key:
      VecGetArray*      -> the Vec argument (arg 0)
      DMDAVecGetArray*  -> the (DM, Vec) argument pair (args 0,1)
  - Matches are tracked as a stack per key per function, in source order,
    so nested/re-entrant Get calls on the same object are handled correctly.
"""

# Map each "Get" function to its matching "Restore" function, and tell the
# script how many of the leading arguments form the identity key
# (1 for Vec* functions, 2 for DMDAVec* functions).
const PAIRS = Dict(
    "VecGetArray"             => ("VecRestoreArray", 1),
    "VecGetArrayRead"         => ("VecRestoreArrayRead", 1),
    "VecGetArrayWrite"        => ("VecRestoreArrayWrite", 1),
    "VecGetArrayPair"         => ("VecRestoreArrayPair", 2),
    "DMDAVecGetArray"         => ("DMDAVecRestoreArray", 2),
    "DMDAVecGetArrayRead"     => ("DMDAVecRestoreArrayRead", 2),
    "DMDAVecGetArrayWrite"    => ("DMDAVecRestoreArrayWrite", 2),
    "DMDAVecGetArrayDOF"      => ("DMDAVecRestoreArrayDOF", 2),
    "DMDAVecGetArrayDOFRead"  => ("DMDAVecRestoreArrayDOFRead", 2),
    "DMDAVecGetArrayDOFWrite" => ("DMDAVecRestoreArrayDOFWrite", 2),
)

const RESTORE_TO_GET = Dict(v[1] => k for (k, v) in PAIRS)
const ALL_NAMES = union(Set(keys(PAIRS)), Set(keys(RESTORE_TO_GET)))
# Longest-first so common prefixes (e.g. VecGetArray vs VecGetArrayRead)
# don't matter for matching order (PCRE backtracks anyway, but this keeps
# it tidy and fast).
const ALL_NAMES_SORTED = sort(collect(ALL_NAMES); by = length, rev = true)

const CALL_RE = Regex("\\b(" * join(ALL_NAMES_SORTED, "|") * ")\\s*\\(")
const COMMENT_RE = r"//.*?$|/\*.*?\*/"sm

function strip_comments(text::AbstractString)::String
    io = IOBuffer()
    last_pos = firstindex(text)
    for m in eachmatch(COMMENT_RE, text)
        write(io, text[last_pos:prevind(text, m.offset)])
        s = m.match
        # Block comments can span multiple lines: replace with the same
        # number of newlines so every later line number stays correct.
        # Line comments never contain a newline, so an empty string is safe.
        if startswith(s, "/*")
            write(io, "\n"^count(==('\n'), s))
        end
        last_pos = m.offset + ncodeunits(m.match)
    end
    write(io, text[last_pos:end])
    return String(take!(io))
end

function split_top_level_args(arg_str::AbstractString)::Vector{String}
    args = String[]
    depth = 0
    cur = IOBuffer()
    for ch in arg_str
        if ch == '('
            depth += 1
            write(cur, ch)
        elseif ch == ')'
            depth -= 1
            write(cur, ch)
        elseif ch == ',' && depth == 0
            push!(args, strip(String(take!(cur))))
            cur = IOBuffer()
        else
            write(cur, ch)
        end
    end
    last_arg = strip(String(take!(cur)))
    if !isempty(last_arg)
        push!(args, last_arg)
    end
    return args
end

function find_matching_paren(text::AbstractString, open_idx::Int)::Int
    depth = 0
    n = ncodeunits(text)
    for i in open_idx:n
        c = text[i]
        if c == '('
            depth += 1
        elseif c == ')'
            depth -= 1
            if depth == 0
                return i
            end
        end
    end
    return -1
end

"""Return a Vector of (start, stop) char offsets of each top-level {...}
block that looks like a function body (heuristic: preceded by ')' before
the '{', at brace depth 0)."""
function iter_functions(text::AbstractString)::Vector{Tuple{Int,Int}}
    result = Tuple{Int,Int}[]
    depth = 0
    start = 0  # 0 means "None"
    n = ncodeunits(text)
    for i in 1:n
        c = text[i]
        if c == '{'
            if depth == 0
                j = i - 1
                while j >= 1 && (text[j] == ' ' || text[j] == '\t' || text[j] == '\r' || text[j] == '\n')
                    j -= 1
                end
                start = (j >= 1 && text[j] == ')') ? i : 0
            end
            depth += 1
        elseif c == '}'
            depth -= 1
            if depth == 0 && start != 0
                push!(result, (start, i))
                start = 0
            end
        end
    end
    return result
end

function key_for(args::Vector{String}, n_key_args::Int)::Vector{String}
    key_args = args[1:min(n_key_args, length(args))]
    return String[String(strip(lstrip(strip(a), '&'))) for a in key_args]
end

function check_function_body(body::AbstractString, filename::AbstractString, base_line::Int, report::Function)
    stack = Dict{Tuple{String,Vector{String}},Vector{Tuple{String,Int}}}()

    for m in eachmatch(CALL_RE, body)
        name = String(m.captures[1])
        open_paren = m.offset + ncodeunits(m.match) - 1  # index of the '(' char
        close_paren = find_matching_paren(body, open_paren)
        close_paren == -1 && continue
        arg_str = body[open_paren+1:close_paren-1]
        args = split_top_level_args(arg_str)
        line = base_line + count(==('\n'), body[firstindex(body):prevind(body, m.offset)])

        if haskey(PAIRS, name)
            _, n_key_args = PAIRS[name]
            prefix = String(split(name, "Get")[1])
            key = (prefix, key_for(args, n_key_args))
            push!(get!(() -> Tuple{String,Int}[], stack, key), (name, line))
        elseif haskey(RESTORE_TO_GET, name)
            get_name = RESTORE_TO_GET[name]
            n_key_args = PAIRS[get_name][2]
            prefix = String(split(get_name, "Get")[1])
            key = (prefix, key_for(args, n_key_args))
            if haskey(stack, key) && !isempty(stack[key])
                pop!(stack[key])
            else
                report(filename, line, "$(name)() with no matching preceding Get call")
            end
        end
    end

    for (_, pending) in stack
        for (get_name, line) in pending
            report(filename, line, "$(get_name)() never matched by a Restore call")
        end
    end
end

function scan_file(path::AbstractString, report::Function)
    local raw
    try
        raw = read(path)
    catch
        return
    end
    # The scanning routines below index the text one byte at a time (like
    # Python's character-indexed strings), which only stays valid if every
    # character is a single byte. Drop any non-ASCII byte (including every
    # byte of a multi-byte UTF-8 sequence) so offsets never land mid-codepoint.
    text = String(UInt8[b for b in raw if b < 0x80])
    text = strip_comments(text)
    for (start, stop) in iter_functions(text)
        base_line = count(==('\n'), text[firstindex(text):prevind(text, start)]) + 1
        check_function_body(SubString(text, start, stop), string(path), base_line, report)
    end
end

function main()
    if length(ARGS) < 1
        println(DOCSTRING)
        exit(1)
    end

    issues = Tuple{String,Int,String}[]
    report(filename, line, msg) = push!(issues, (String(filename), line, String(msg)))

    exts = Set([".c", ".h", ".cxx", ".cpp", ".hpp", ".cu", ".cuh"])

    for root in ARGS
        if isfile(root)
            scan_file(root, report)
        elseif isdir(root)
            for (dirpath, _, files) in walkdir(root)
                for f in files
                    if splitext(f)[2] in exts
                        scan_file(joinpath(dirpath, f), report)
                    end
                end
            end
        end
    end

    sort!(issues)
    for (filename, line, msg) in issues
        println("$(filename):$(line): $(msg)")
    end

    println()
    println("$(length(issues)) potential mismatch(es) found.")
    exit(isempty(issues) ? 0 : 1)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
