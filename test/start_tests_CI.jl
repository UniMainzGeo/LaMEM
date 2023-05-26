# this starts the testing suite from the CI
# the difference is that we use LaMEM compiled vs dynamic libraries so we need to add the environment
using Pkg

# run test suite
Pkg.test("LaMEM_C", test_args=["use_dynamic_lib"])


exit()
