#
#   Copyright 2022-2023 LAIKA. Authored by Mitch Prater.
# 
#   Licensed under the Apache License Version 2.0 http://apache.org/licenses/LICENSE-2.0,
#   or the MIT license http://opensource.org/licenses/MIT, at your option.
# 
#   This program may not be copied, modified, or distributed except according to those terms.
#

#
# PIXAR_ROOT and RMAN_VERSION are required to be set.
#
ifndef PIXAR_ROOT
    $(error PIXAR_ROOT has not been set. This is required for the build system to function.)
endif
ifndef RMAN_VERSION
    $(error RMAN_VERSION has not been set. This must be done by the top-level Makefile, the environment, or in the 'make' command)
endif

# Set RMANTREE is if isn't already.
RMANTREE ?= $(PIXAR_ROOT)/RenderManProServer-$(RMAN_VERSION)


#------------------------------------------------------------------
# Make functionality.
#	.so files are made in-place from their .cpp file.
#	The .so files are then installed to $(install_dir)
#	with their flattened base name for USD compatibility.
#------------------------------------------------------------------
#
# SUBDIRS and EXCLUDE control variables determine which .cpp files are made.
#
# These directories will be searched for .cpp files.
SUBDIRS := bxdf

# These shaders won't be built: directory/basename
EXCLUDE := 

# SUBDIRS and EXCLUDE functionality.
#	$(src) is the list of .cpp files to build from.
#	$(dso) is the list of .so files to build.
src := $(foreach dir,$(SUBDIRS),$(wildcard $(dir)/*.cpp))

EXCLUDE := $(basename $(EXCLUDE))
EXCLUDE := $(addsuffix .cpp, $(EXCLUDE))

src := $(filter-out $(EXCLUDE), $(src))
args := $(foreach file,$(src),$(addprefix $(dir $(file)),Args/$(notdir $(file:%.cpp=%.args))))

obj := $(src:%.cpp=%.o)
dso := $(src:%.cpp=%.so)

# Relevant make system directories.
SRCDIR ?= ..
PYTHONDIR ?= $(SRCDIR)/python3

# Location to install the flattened shaders:
# $(SRCDIR)/build/
# └── $(RMAN_VERSION)
#     └── plugins
#         ├── flattened_ShaderName.so
#         ├── ...
#         └── Args
#             ├── flattened_ShaderName.args
#             └── ...
DSTDIR ?= $(SRCDIR)/build
install_dir := $(DSTDIR)/${RMAN_VERSION}/plugins

#
# Default goal: the target of 'make'.
#
.DEFAULT_GOAL = all

# Targets whose time stamps we want to ignore.
.PHONY : all clean install

# The default target's prerequisites.
all : install $(dso) $(args)

# The installation requires $(install_dir)/Args to exist.
install : $(install_dir)/Args

# Don't worry if $(install_dir) already exists, just make it and Args.
$(install_dir)/Args : 
	@ mkdir -p $(install_dir)/Args

#
# Compiler settings.
#
CLANG = clang++

CC = $(CLANG) -fdiagnostics-color -ferror-limit=5
LD = $(CLANG) -fdiagnostics-color -shared

CFLAGS = -Wall -Wno-switch -Wno-reorder -Wno-unused-private-field -fPIC -std=c++17 -O3
LFLAGS = 

INCLUDE = -I. -I./include -I${RMANTREE}/include -I${RMANTREE}/plugins

#
# Plugin build rules.
#
# include directory header file dependencies.
includeFiles := ./include/*.h
$(so): $(includeFiles)

# Individual dependencies.
bxdf/PraterFuzz.o: bxdf/PraterFuzzSampling.inl
bxdf/PraterScatter.o: bxdf/PraterScatterSampling.inl
bxdf/KPVelvet.o: bxdf/KPVelvetSampling.inl
bxdf/HGScatter.o: bxdf/HGScatterSampling.inl

# Make the .so file in-place and copy it to the
# $(install_dir) using its flattened file name.
%.so : %.o Args/%.args
	$(LD) $(LFLAGS) $(INCLUDE) -o $@ $<
	env python3 $(PYTHONDIR)/install_shaders.py --copy $(install_dir) $@

%.o : %.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ -c $<

#
# Clean rules.
#
clean_local:
	@ echo "cpp directory clean."
	@ -rm -f $(shell find -L . -name "*.o" -print)
	@ -rm -f $(shell find -L . -name "*.so" -print)

clean: clean_local
	@ -rm -f $(shell find -L $(install_dir) -name "*.so" -print)
	@ -rm -f $(shell find -L $(install_dir) -name "*.args" -print)


#
# Helpful rules.
#
help:
	@ echo "------------------------------------------------------------------------"
	@ echo "PIXAR_ROOT must be set to the location of the RenderMan installation:"
	@ echo "e.g. /opt/pixar"
	@ echo "RMAN_VERSION must be set to the RenderMan version you wish to make:"
	@ echo "e.g. 25.2"
	@ echo ""
	@ echo "Current settings:"
	@ echo "PIXAR_ROOT: $(PIXAR_ROOT)"
	@ echo "RMAN_VERSION: $(RMAN_VERSION)"
	@ echo "RMANTREE:   $(RMANTREE)"
	@ echo "CURDIR:     $(CURDIR)"
	@ echo "PYTHONDIR:  $(PYTHONDIR)"
	@ echo "DSTDIR:     $(DSTDIR)"
	@ echo "build dir:  $(install_dir)"
	@ echo "------------------------------------------------------------------------"

show_src:
	@ echo "The source .cpp files to compile:"
	@ echo $(src) | sed "s/ /\\n/g" | sort
	@ echo ""

show_args:
	@ echo "The .args files:"
	@ echo $(args) | sed "s/ /\\n/g" | sort
	@ echo ""

show_so:
	@ echo "The target .so files to build:"
	@ echo $(dso) | sed "s/ /\\n/g" | sort
	@ echo ""
	@ echo "The built .so files that exist in the ris directory tree:"
	@ find -L . -name "*.so" -print | sort
	@ echo ""

show_build:
	@ echo "The installed .so files:"
	@ find -L $(install_dir) -name "*.so" -print | sort
	@ echo ""
