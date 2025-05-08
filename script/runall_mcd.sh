#!/usr/bin/env bash
catkin_make -C /home/$USER/work/ctemlo_WS;
source /home/$USER/work/ctemlo_WS/devel/setup.bash # TODO: Please change this one based on your env
export EPOC_DIR=/mnt/e/bag/MCD/CTE-MLO # Log path
export DATASET_LOCATION=/mnt/e/bag/MCD # bag path


# Get the current directory
CURR_DIR=$(pwd)
# Get the location of the viral package
roscd cte_mlo
PACKAGE_DIR=$(pwd)
# Return to the current dir, print the directions
cd $CURR_DIR
echo CURRENT DIR: $CURR_DIR
echo VIRAL DIR:   $PACKAGE_DIR
export CAPTURE_SCREEN=false;
export LOG_DATA=true;

#region 0 UWB NO VIS --------------------------------------------------------------------------------------------------

# export EPOC_DIR=$1;
# export DATASET_LOCATION=$2;
# export ROS_PKG_DIR=$3;
# export EXP_NAME=$4;
# export CAPTURE_SCREEN=$5;
# export LOG_DATA=$6;
# export LOG_DUR=$7;
# export FUSE_UWB=$8;
# export FUSE_VIS=$9;
# export UWB_BIAS=${10};
# export ANC_ID_MAX=${11};

wait;
./run_one_bag_mcd.sh $EPOC_DIR $DATASET_LOCATION $PACKAGE_DIR ntu_day_01 $CAPTURE_SCREEN $LOG_DATA 450 0 1 0.75 -1;
wait;
./run_one_bag_mcd.sh $EPOC_DIR $DATASET_LOCATION $PACKAGE_DIR ntu_day_02 $CAPTURE_SCREEN $LOG_DATA 450 0 1 0.75 -1;
wait;
./run_one_bag_mcd.sh $EPOC_DIR $DATASET_LOCATION $PACKAGE_DIR ntu_day_10 $CAPTURE_SCREEN $LOG_DATA 450 0 1 0.75 -1;

wait;
./run_one_bag_mcd.sh $EPOC_DIR $DATASET_LOCATION $PACKAGE_DIR ntu_night_04 $CAPTURE_SCREEN $LOG_DATA 450 0 1 0.75 -1;
wait;
./run_one_bag_mcd.sh $EPOC_DIR $DATASET_LOCATION $PACKAGE_DIR ntu_night_08 $CAPTURE_SCREEN $LOG_DATA 450 0 1 0.75 -1;
wait;
./run_one_bag_mcd.sh $EPOC_DIR $DATASET_LOCATION $PACKAGE_DIR ntu_night_13 $CAPTURE_SCREEN $LOG_DATA 450 0 1 0.75 -1;

#endregion NO UWB NO VIS ----------------------------------------------------------------------------------------------


#region ## Poweroff ---------------------------------------------------------------------------------------------------

# wait;
# poweroff

#endregion ## Poweroff ------------------------------------------------------------------------------------------------
