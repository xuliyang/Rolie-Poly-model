#!/bin/sh

#. $WM_PROJECT_DIR/bin/tools/CleanFunctions

#cleanCase



#PROGRAM="./xly_Oldroyd-B"
PROGRAM="Rolie-Poly-CyM"
PNUM="12"


    TIME=`date +%m%d%k%M`
    RUN="mpiexec -hosts $1 -n $PNUM $PROGRAM -parallel"
    echo $RUN
    $RUN 2>&1 > run-${TIME}.log

    echo ""
    echo "Done"
