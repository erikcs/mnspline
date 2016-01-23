# compile to a matlab executable on OSX with OpenMP support

# System dependent parameters that may have to be changed:
CC='gcc-5'
MATLABDIR='/Applications/MATLAB_R2015b.app'
XCODEINC='/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.11.sdk'
OMPI='/usr/local/include/libiomp'
OMPL='/usr/local/lib'
OSXV=10.11

#############################################################################

SOURCE='mnspline_mex'
DEPS='mnspline'
TARGET='mnspline.mexmaci64'

OPTFLAGS='-O3'

$CC -c -DMX_COMPAT_32   -DMATLAB_MEX_FILE -I$OMPI \
        -I"$MATLABDIR/extern/include" -I"$MATLABDIR/simulink/include" \
        \-fno-common -arch x86_64 -mmacosx-version-min=$OSXV\
         -fexceptions -isysroot $XCODEINC \
        -fopenmp $OPTFLAGS -DNDEBUG "$SOURCE.c" "../src/$DEPS.c" 

$CC \-Wl,-twolevel_namespace -undefined error -arch x86_64\
        -mmacosx-version-min=$OSXV -Wl,-syslibroot,$XCODEINC -bundle\
        -Wl,-exported_symbols_list,"$MATLABDIR/extern/lib/maci64/mexFunction.map"\
        -fopenmp "$SOURCE.o" "$DEPS.o" -O -Wl,-exported_symbols_list,\
        "$MATLABDIR/extern/lib/maci64/mexFunction.map"\
        -L$OMPL   -L"$MATLABDIR/bin/maci64" -lmx -lmex -lmat -lc++\
        -o $TARGET

rm "$SOURCE.o"
rm "$DEPS.o"
