AC_DEFUN([ACX_CHECK_CXX_FLAGS],
[
AC_REQUIRE([AC_PROG_CXX])
AC_CACHE_CHECK(whether ${CXX-cc} accepts $1, ac_cv_$2,
[echo 'void f(){}' > conftest.c
if test -z "`${CXX-cc} $1 -c conftest.c 2>&1`"; then
	ac_cv_$2=yes
else
	ac_cv_$2=no
fi
rm -f conftest*
])
if test "$ac_cv_$2" = yes; then
	:
	$3
else
	:
	$4
fi
])

AC_DEFUN([ACX_PROG_GXX_VERSION],
[
AC_REQUIRE([AC_PROG_CXX])
AC_CACHE_CHECK(whether we are using g++ $1.$2 or later, ac_cv_prog_gxx_$1_$2,
[
dnl The semicolon after "yes" below is to pacify NeXT's syntax-checking cpp.
cat > conftest.c <<EOF
#ifdef __GNUG__
#  if (__GNUC__ > $1) || (__GNUC__ == $1 && __GNUC_MINOR__ >= $2)
     yes;
#  endif
#endif
EOF
if AC_TRY_COMMAND(${CXX-cc} -E conftest.c) | egrep yes >/dev/null 2>&1; then
  ac_cv_prog_gxx_$1_$2=yes
else
  ac_cv_prog_gxx_$1_$2=no
fi
])
if test "$ac_cv_prog_gxx_$1_$2" = yes; then
	:
	$3
else
	:
	$4
fi
])

dnl AC_DEFUN(ACX_PROG_CC_MAXOPT,
dnl [
dnl AC_REQUIRE([AC_PROG_CC])
dnl AC_REQUIRE([ACX_PROG_CC_EGCS])
dnl AC_REQUIRE([AC_CANONICAL_HOST])
dnl 
dnl # Try to determine "good" native compiler flags if none specified on command
dnl # line
dnl if test "$ac_test_CFLAGS" != "set"; then
dnl   CFLAGS=""
dnl   case "${host_cpu}-${host_os}" in
dnl 
dnl   *linux*)
dnl 	echo "*******************************************************"
dnl 	echo "*       Congratulations! You are running linux.       *"
dnl 	echo "*******************************************************"
dnl 	;;
dnl   sparc-solaris2*) if test "$CC" = cc; then
dnl                     CFLAGS="-native -fast -xO5 -dalign"
dnl                  fi;;
dnl 
dnl   alpha*-osf*)  if test "$CC" = cc; then
dnl                     CFLAGS="-newc -w0 -O5 -ansi_alias -ansi_args -fp_reorder -tune host -arch host -std1"
dnl                 fi;;
dnl 
dnl   hppa*-hpux*)  if test "$CC" = cc; then
dnl                     CFLAGS="-Ae +O3 +Oall"
dnl                 fi;;
dnl 
dnl    *-aix*)
dnl 	if test "$CC" = cc -o "$CC" = xlc; then
dnl 		ACX_CHECK_CC_FLAGS([-qarch=auto -qtune=auto], qarch_auto,
dnl 			[CFLAGS="-O3 -qansialias -w -qarch=auto -qtune=auto"],
dnl 			[CFLAGS="-O3 -qansialias -w"
dnl 	echo "*******************************************************"
dnl 	echo "*  You seem to have AIX and the IBM compiler.  It is  *"
dnl 	echo "*  recommended for best performance that you use:     *"
dnl 	echo "*                                                     *"
dnl 	echo "*    CFLAGS=-O3 -qarch=xxx -qtune=xxx -qansialias -w  *"
dnl 	echo "*                      ^^^        ^^^                 *"
dnl 	echo "*  where xxx is pwr2, pwr3, 604, or whatever kind of  *"
dnl         echo "*  CPU you have.  (Set the CFLAGS environment var.    *"
dnl         echo "*  and re-run configure.)  For more info, man cc.     *"
dnl 	echo "*******************************************************"
dnl 			])
dnl         fi;;
dnl   esac
dnl 
dnl   # use default flags for gcc on all systems
dnl   if test $ac_cv_prog_gcc = yes; then
dnl      CFLAGS="-O3 -fomit-frame-pointer"
dnl   fi
dnl 
dnl   # the egcs scheduler is too smart and destroys our own schedule.
dnl   # Disable the first instruction scheduling pass.  The second
dnl   # scheduling pass (after register reload) is ok.
dnl   if test "$acx_prog_egcs" = yes; then
dnl      CFLAGS="$CFLAGS -fno-schedule-insns -fschedule-insns2"
dnl   fi
dnl 
dnl   # test for gcc-specific flags:
dnl   if test $ac_cv_prog_gcc = yes; then
dnl     # -malign-double for x86 systems
dnl     ACX_CHECK_CC_FLAGS(-malign-double,align_double,
dnl 	CFLAGS="$CFLAGS -malign-double")
dnl     # -fstrict-aliasing for gcc-2.95+
dnl     ACX_CHECK_CC_FLAGS(-fstrict-aliasing,fstrict_aliasing,
dnl 	CFLAGS="$CFLAGS -fstrict-aliasing")
dnl   fi
dnl 
dnl   CPU_FLAGS=""
dnl   if test "$GCC" = "yes"; then
dnl 	  dnl try to guess correct CPU flags, at least for linux
dnl 	  case "${host_cpu}" in
dnl 	  i586*)  ACX_CHECK_CC_FLAGS(-mcpu=pentium,cpu_pentium,
dnl 			[CPU_FLAGS=-mcpu=pentium],
dnl 			[ACX_CHECK_CC_FLAGS(-mpentium,pentium,
dnl 				[CPU_FLAGS=-mpentium])])
dnl 		  ;;
dnl 	  i686*)  ACX_CHECK_CC_FLAGS(-mcpu=pentiumpro,cpu_pentiumpro,
dnl 			[CPU_FLAGS=-mcpu=pentiumpro],
dnl 			[ACX_CHECK_CC_FLAGS(-mpentiumpro,pentiumpro,
dnl 				[CPU_FLAGS=-mpentiumpro])])
dnl 		  ;;
dnl 	  alphaev4-*)  ACX_CHECK_CC_FLAGS(-mcpu=ev4,cpu_ev4,
dnl 			[CPU_FLAGS=-mcpu=ev4])
dnl 		  ;;
dnl 	  alphaev56-*)  ACX_CHECK_CC_FLAGS(-mcpu=ev56,cpu_ev56,
dnl 			[CPU_FLAGS=-mcpu=ev56])
dnl 		  ;;
dnl 	  alphaev5-*)  ACX_CHECK_CC_FLAGS(-mcpu=ev5,cpu_ev5,
dnl 			[CPU_FLAGS=-mcpu=ev5])
dnl 		  ;;
dnl 	  alphaev6-*)  ACX_CHECK_CC_FLAGS(-mcpu=ev6,cpu_ev6,
dnl 			[CPU_FLAGS=-mcpu=ev6])
dnl 		  ;;
dnl 	  powerpc*)
dnl 		cputype=`(grep cpu /proc/cpuinfo | head -1 | cut -d: -f2 | sed 's/ //g') 2> /dev/null`
dnl 		is60x=`echo $cputype | egrep "^60[0-9]e?$"`
dnl 		if test -n "$is60x"; then
dnl 			ACX_CHECK_CC_FLAGS(-mcpu=$cputype,m_cpu_60x,
dnl 				CPU_FLAGS=-mcpu=$cputype)
dnl 		elif test "$cputype" = 750; then
dnl                         ACX_PROG_GCC_VERSION(2,95,
dnl                                 ACX_CHECK_CC_FLAGS(-mcpu=750,m_cpu_750,
dnl 					CPU_FLAGS=-mcpu=750))
dnl 		fi
dnl 		if test -z "$CPU_FLAGS"; then
dnl 		        ACX_CHECK_CC_FLAGS(-mcpu=powerpc,m_cpu_powerpc,
dnl 				CPU_FLAGS=-mcpu=powerpc)
dnl 		fi
dnl 		if test -z "$CPU_FLAGS"; then
dnl 			ACX_CHECK_CC_FLAGS(-mpowerpc,m_powerpc,
dnl 				CPU_FLAGS=-mpowerpc)
dnl 		fi
dnl 	  esac
dnl   fi
dnl 
dnl   if test -n "$CPU_FLAGS"; then
dnl         CFLAGS="$CFLAGS $CPU_FLAGS"
dnl   fi
dnl 
dnl   if test -z "$CFLAGS"; then
dnl 	echo ""
dnl 	echo "********************************************************"
dnl         echo "* WARNING: Don't know the best CFLAGS for this system  *"
dnl         echo "* Use  make CFLAGS=..., or edit the top level Makefile *"
dnl 	echo "* (otherwise, a default of CFLAGS=-O3 will be used)    *"
dnl 	echo "********************************************************"
dnl 	echo ""
dnl         CFLAGS="-O3"
dnl   fi
dnl 
dnl   ACX_CHECK_CC_FLAGS(${CFLAGS}, guessed_cflags, , [
dnl 	echo ""
dnl         echo "********************************************************"
dnl         echo "* WARNING: The guessed CFLAGS don't seem to work with  *"
dnl         echo "* your compiler.                                       *"
dnl         echo "* Use  make CFLAGS=..., or edit the top level Makefile *"
dnl         echo "********************************************************"
dnl         echo ""
dnl         CFLAGS=""
dnl   ])
dnl 
dnl fi
dnl ])
