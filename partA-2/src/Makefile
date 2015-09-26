# The most recent VXL library is available in the following path:
VXLDIR = /h/u1/mangas/public/vxl

# Include directories and flags for compilation
CFLAGS = -I$(VXLDIR)/ -I$(VXLDIR)/vxl-1.17.0 -I$(VXLDIR)/vxl-1.17.0/core -I$(VXLDIR)/vxl-1.17.0/vcl -I$(VXLDIR)/bin/vcl -I$(VXLDIR)/bin/core -I$(VXLDIR)/vxl-1.17.0/core/vil -Wno-deprecated -DJPEG_LIB_VERSION=80

CC = g++

#fltk library paths and libraries
FLTKLDFLAGS = -L/local/lib -L/usr/lib/x86_64-linux-gnu/ -Wl,-rpath,/local/lib -lfltk_gl -lGLU -lGL -lfltk -L/usr/X11R6/lib -lpthread -lm -lXext -lX11 -lsupc++

#vxl library paths and libraries
VXLLDFLAGS = -L$(VXLDIR)/bin/lib -lvil -lvnl_algo -lvnl_io -lvnl -lv3p_netlib -lvcl -lvul -lpng -ltiff -lopenjpeg2 -lgeotiff
LDFLAGS = -lstdc++ $(FLTKLDFLAGS) $(VXLLDFLAGS)

FLUID = fluid

%.o: %.cxx
	$(CC) -c $(CFLAGS) $(CPPFLAGS) -o $@ $<

%.cxx: %.fl
	$(FLUID) -c $<

BASIC_OBJ   = gl/glutils.o gl/Texture.o main.o file/load_image.o

IMDRAW_OBJ  = imdraw/imdraw_utils.o imdraw/imdraw_init.o imdraw/imdraw_draw.o imdraw/imdraw_handle.o imdraw/read_drawing.o

UI_OBJ      = VisCompUI.o

UI_CPP      = VisCompUI.cxx VisCompUI.h

MATTING_OBJ = matting/matting.o

STUDENT_OBJ = matting/matting_algorithm.o

imdraw-code/viscomp:	$(UI_CPP) $(BASIC_OBJ) $(UI_OBJ) $(STUDENT_OBJ) $(MATTING_OBJ) $(IMDRAW_OBJ)

	$(CC) -o ../bin/viscomp $(BASIC_OBJ) $(UI_OBJ) $(MATTING_OBJ) $(IMDRAW_OBJ) $(STUDENT_OBJ) $(LDFLAGS)

clean:		

	rm -rf $(BASIC_OBJ) $(UI_OBJ) $(STUDENT_OBJ) $(UI_CPP) $(MATTING_OBJ) $(IMDRAW_OBJ) 

