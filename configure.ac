AC_PREREQ([2.65])

AC_INIT("CNEr")

## Checks for libraries


## Checks for headers
AC_CHECK_HEADERS([arpa/inet.h netinet/in.h signal.h stdio.h stdlib.h 
                  stdarg.h string.h ctype.h sys/types.h sys/stat.h 
                  strings.h fcntl.h assert.h setjmp.h sys/time.h math.h 
                  errno.h unistd.h libgen.h pthread.h sys/socket.h 
                  netdb.h utime.h limits.h sys/utsname.h pwd.h termios.h 
                  sys/wait.h regex.h dirent.h sys/file.h])

## Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_OFF_T
AC_TYPE_PID_T
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T
AC_TYPE_UID_T

## Checks for library functions.
AC_FUNC_FORK
AC_FUNC_MALLOC
AC_FUNC_MKTIME
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([bzero dup2 getcwd gethostbyname gettimeofday memmove 
                memset mkdir regcomp select socket sqrt strcasecmp strchr 
                strcspn strdup strerror strrchr strspn strstr uname])

AC_OUTPUT(src/Makevars)
