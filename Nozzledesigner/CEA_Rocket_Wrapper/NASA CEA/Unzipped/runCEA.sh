#!/bin/sh
# 02/12/2004  Minna.M.Chao@nasa.gov 216-433-5197
# J2SDK 1.4.1_02 and Visual Cafe Enterprise 4.52
# %% %%
# CEA_HOME=/u/CEAgui
curDir=`pwd`
which java
java -version

# JDK_HOME = /usr/jdk1.4/; export JDK_HOME
# PATH = $JDK_HOME/bin:$PATH; export PATH
CLASSPATH="CEAgui.jar"

echo "PWD=" `pwd`
echo " "
echo "CLASSPATH=" $CLASSPATH

echo " "
echo "PATH=" $PATH

echo " "
java -version

echo " "
CMD="java  -classpath CEAgui.jar CEAgui"
#CMD=java  -classpath ${CLASSPATH} CEAgui

echo ${CMD}
${CMD}
