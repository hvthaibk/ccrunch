set(CTEST_PROJECT_NAME          "ccrunch")
set(CTEST_NIGHTLY_START_TIME    "01:00:00 UTC")

set(CTEST_DROP_METHOD           "http")
set(CTEST_DROP_SITE             "my.cdash.org")
set(CTEST_DROP_LOCATION         "/submit.php?project=ccrunch")
set(CTEST_DROP_SITE_CDASH       TRUE)

set(CTEST_SOURCE_DIRECTORY      "${PROJECT_SOURCE_DIR}")
set(CTEST_BINARY_DIRECTORY      "${PROJECT_SOURCE_DIR}/build64")

find_program(HOSTNAME_CMD NAMES hostname)
exec_program(${HOSTNAME_CMD} ARGS OUTPUT_VARIABLE HOSTNAME)

find_program(UNAME NAMES uname)
macro(getuname name flag)
  exec_program("${UNAME}" ARGS "${flag}" OUTPUT_VARIABLE "${name}")
endmacro(getuname)

getuname(osname -s)
getuname(osrel  -r)
getuname(cpu    -m)

set(SITE            "${HOSTNAME}")
set(BUILDNAME       "${osname}--${osrel}--${cpu}--GCC${GCC_VERSION}")

set(WITH_MEMCHECK TRUE)
set(WITH_COVERAGE TRUE)

set(CTEST_COMMAND
   "/usr/bin/ctest -D Experimental"
   )
#   "/usr/bin/ctest -D NightlyStart -D NightlyUpdate -D NightlyConfigure -D NightlyBuild -D NightlyTest -D NightlySubmit"
#   "/usr/bin/ctest -D NightlyMemCheck -D NightlySubmit"
