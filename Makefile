
#Compiler
CXX = g++

#Debug or Release
DEBUG = -g
RELEASE = -O2
EXEC_TYPE = ${DEBUG}

#Compiler options
CXXFLAGS = -Wall ${EXEC_TYPE} $(shell root-config --cflags)

#PATH
#PATH_X = -I/usr/X11R6/lib
#PATH_OPENGL = -I/usr/include/GL
#PATH_BOOST = -I /usr/include/boost
PATH_ITS = -I./includes/
INCPATH = ${PATH_ITS}

#FLAGS
#FLAGS_X = -lX11 -lXmu -lXi 
#FLAGS_OPENGL = -lGL -lGLU
#FLAGS_BOOST = -lboost_system -lboost_filesystem
OTHER_FLAGS = -lm -lstdc++ -std=c++11
LDFLAGS =  ${OTHER_FLAGS} $(shell root-config --glibs) 

#executable name
EXEC = tb2017ana

SOURCE_DIR = source

#list here all the source files
SRC = main.cpp $(SOURCE_DIR)/histos.cpp $(SOURCE_DIR)/MAMBA_presetting.cpp $(SOURCE_DIR)/STRAW_presetting.cpp $(SOURCE_DIR)/linking_data.cpp 

#A directory for build artefacts
BUILD_DIR = Build

#OBJ files (if sources in cpp)
OBJ = $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(SRC))

#what we are trying to build
all: $(EXEC)

#linkage
$(EXEC): ${OBJ}
	@echo 
	@echo =============================== [LINKAGE ...] =============================== 
	@echo

	$(CC) -o $(BUILD_DIR)/$@ $^ $(LDFLAGS)
	@echo
	@echo ======================= [progam compiled succesfully] =======================
	@echo 

$(OBJ): | ${BUILD_DIR}


$(BUILD_DIR):
	@mkdir -p $@

#compile every source files
$(BUILD_DIR)/%.o : %.cpp

	@echo --------- build [$<] ----------

	@mkdir -p $(@D) # http://www.gnu.org/software/make/manual/make.html#index-_0024_0028_0040D_0029

	$(CXX) $(CXXFLAGS) $(INCPATH) -c $<  -o $@
	@echo OK [$<] - [$@]
	@echo



#make clean
clean:
	rm -rf $(BUILD_DIR)
