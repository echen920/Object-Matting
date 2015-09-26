CC = gcc
CFLAGS = `fltk-config --cxxflags`
#LDFLAGS = `fltk-config --use-gl --ldstaticflags` -lstdc++
LDFLAGS = -L/local/lib -Wl,-rpath,/local/lib -lfltk_gl -lGLU -lGL -lfltk -L/usr/X11R6/lib -lpthread -lm -lXext -lX11 -lsupc++

%.o: %.cxx
	$(CC) -c $(CFLAGS) $(CPPFLAGS) -o $@ $<

CubeView: CubeView.o CubeViewUI.o CubeMain.o
	$(CC) -o CubeView CubeView.o CubeViewUI.o CubeMain.o $(LDFLAGS)

hello: hello.o
	$(CC) -o hello hello.o $(LDFLAGS)
