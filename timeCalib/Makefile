CXX=`root-config --cxx`
CXXFLAGS=`root-config --cflags`
LDFLAGS=`root-config --ldflags`
LDLIBS=`root-config --glibs`

all: timeCalib

% : %.C
	g++ -g -Wall `root-config --cflags --libs` -L$(ROOTSYS)/lib $? -o $@

clean:
	rm -f *.o
	rm -f $(TARGET)





