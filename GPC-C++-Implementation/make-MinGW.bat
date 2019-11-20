@ECHO OFF
MKDIR build
CD build
cmake -G "Unix Makefiles" ..
make

