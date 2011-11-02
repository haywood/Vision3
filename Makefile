# Makefile for Computer Vision HW5
# Comments begin with a # (like this line)

# Your UNI
UNI = abc1234

# VARIOUS VARIABLES (you will need to change some of these)
# The compiler. If you're writing code using c++, change gcc to g++
CC = g++

# Compile time flags. If you want to compile with debugging turned off, remove the -g
CFLAGS = -g

# Libraries to link against. You should only need to link against the math library (-lm).
LFLAGS = -lm

# Enter the names of your source files here. If you have header files as well,
# then do not include them here (they will be autoincluded).
P1SRC = p1.cpp
P2SRC = p2.cpp
P3SRC = p3.cpp
P4SRC = p4.cpp
P5SRC = p5.cpp

# This contains a list of various images and other files generated from 
# running your programs. By default, it contains the filenames of all files
# generated from running the 'test' target. If you generate any other files,
# make sure this is correct. This list of files is deleted when 'make clean'
# or 'make submit' are run (to prevent submitting all these files).
OUTPUTS_TO_CLEAN = sphere.txt directions.txt sphere_mask.pgm vase_mask.pgm sphere_normal.pbm vase_normal.pbm sphere_albedo.pgm vase_albedo.pgm sphere_gradient.dat vase_gradient.dat

# THE BUILD RULES (you shouldn't need to edit any of these)
# 'make test' builds all 5 programs and then tests them.
# You should check to see that the final output images are correct
test: test1 test2 test3 test4 test5
	
# 'make test1' tests only program 1
test1: p1
	./p1 sphere0.pgm sphere.txt

# 'make test2' tests only program 2
test2: p2
	./p2 sphere.txt sphere1.pgm sphere2.pgm sphere3.pgm sphere4.pgm sphere5.pgm directions.txt

# 'make test3' tests only program 3
test3: p3
	./p3 sphere1.pgm sphere2.pgm sphere3.pgm sphere4.pgm sphere5.pgm sphere_mask.pgm
	./p3 vase1.pgm vase2.pgm vase3.pgm vase4.pgm vase5.pgm vase_mask.pgm

# 'make test4' tests only program 4
test4: p4
	./p4 directions.txt sphere1.pgm sphere2.pgm sphere3.pgm sphere4.pgm sphere5.pgm sphere_mask.pgm sphere_normal.ppm sphere_albedo.pgm
	./p4 directions.txt vase1.pgm vase2.pgm vase3.pgm vase4.pgm vase5.pgm vase_mask.pgm vase_normal.ppm vase_albedo.pgm

# 'make test5' tests only program 5
test5: p5
	./p5 sphere_gradient.dat sphere_seed.txt sphere_mask.pgm sphere_depth.pgm
	./p5 vase_gradient.dat vase_seed.txt vase_mask.pgm vase_depth.pgm

# 'make all' builds all 4 programs (but doesn't run them)
all: p1 p2 p3 p4 p5

# The utilities file. Leave this as is.
UTILS = vision_utilities.c
UTILSH = vision_utilities.h

# 'make p1' builds p1 only
p1: $(P1SRC) $(UTILS)
	$(CC) $(CFLAGS) -o $@ $(P1SRC) $(UTILS) $(LFLAGS)

# 'make p2' builds p2 only
p2: $(P2SRC) $(UTILS)
	$(CC) $(CFLAGS) -o $@ $(P2SRC) $(UTILS) $(LFLAGS)

# 'make p3' builds p3 only
p3: $(P3SRC) $(UTILS)
	$(CC) $(CFLAGS) -o $@ $(P3SRC) $(UTILS) $(LFLAGS)

# 'make p4' builds p4 only
p4: $(P4SRC) $(UTILS)
	$(CC) $(CFLAGS) -o $@ $(P4SRC) $(UTILS) $(LFLAGS)

# 'make p5' builds p5 only
p5: $(P5SRC) $(UTILS)
	$(CC) $(CFLAGS) -o $@ $(P5SRC) $(UTILS) $(LFLAGS)

# the rule for building the vision utilities
vision_utilities.o: vision_utilities.c vision_utilities.h
	$(CC) -c -Wall -pedantic -o $@ vision_utilities.c

# 'make submit' creates a tar of all files in the directory
# NOTE 1: It also runs 'make clean' before submission.
 submit: clean
	tar -czvf $(UNI)-hw4.tar.gz $(P1SRC) $(P2SRC) $(P3SRC) $(P4SRC) $(P5SRC) $(UTILS) $(UTILSH) vase_seed.txt sphere_seed.txt Makefile README

# 'make clean' removes all output and temporary files
clean:
	-rm -f p1 p2 p3 p4 p5 vision_utilities.o $(OUTPUTS_TO_CLEAN) core CORE
