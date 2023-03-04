# Inverse Kinematics

This is Assignment 2 of the course DD2410 Introduction to Robotics.

## How to

To run this repository, ubuntu 16.04 with ROS-kinetic is needed.

Download this repository, and put the folder into the ROS workspace, catkin_make and source the devel/setup.bash.

In this repository, I only modify `kinematics_assignment_metapackage/kinematics_assignment/scripts/IK_functions.py` files, and this is also where the algorithm is.

PART 1: SCARA ROBOT
	
To visualize the robot, run: 
	
          roslaunch kinematics_assignment scara_launch.launch
	
PART 2: KUKA ROBOT
	
To visualize the robot, run: 

           roslaunch kinematics_assignment kuka_launch.launch

