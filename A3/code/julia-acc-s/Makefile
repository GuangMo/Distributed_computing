CC=pgcc
CFLAGS=-acc -Minfo=accel -O2 -ta=nvidia,cc13 
LDFLAGS=-acc -ta=nvidia,cc13

OBJS =  main.o julia.o savebmp.o color.o getparams.o 

julia:  $(OBJS)
	$(CC) $(LDFLAGS) -o julia_acc_s $?


# this runs are on Mac. On Linux, e.g. penguin, replace open by gthumb
run0: 
	./julia_acc_s 0 -0.4 0.6  -0.2 -0.1 1  1.1 1000 1000 1000  image.bmp ; photoqt image.bmp &

run1: julia_acc_s
	./julia_acc_s 0 -0.4 0.6  -0.2 -0.1 1  1.1 1000 1000 3000  image.bmp ; photoqt image.bmp &

run2: julia_acc_s
	./julia_acc_s 0 -0.4 0.6  -0.181862 -0.181772 1.019085 1.019175 1000 1000 100  image.bmp ; photoqt image.bmp &

run3: julia_acc_s
	./julia_acc_s 0 -0.4 0.6  -0.181862 -0.181772 1.019085 1.019175 2000 2000 10000  image.bmp ; photoqt image.bmp &

run4: julia_acc_s
	./julia_acc_s 1 -0.8  0.156  -2 2 -2 2 1000 1000 100  image.bmp ; photoqt image.bmp &

run5: julia_acc_s
	./julia_acc_s 1 -0.8  0.156  -2 2 -2 2 1000 1000 1000  image.bmp ; photoqt image.bmp &

run6: julia_acc_s
	./julia_acc_s 1 -0.8  0.156  -.2 .2 -.2 .2 4000 4000 10000  image.bmp ; photoqt image.bmp&

run7: julia_acc_s
	./julia_acc_s 1 0.285 0.01  -2 2 -2 2 1000 1000 1000  image.bmp ; photoqt image.bmp &

run8: julia_acc_s
	./julia_acc_s 1 0.285 0.01  -.2 .2 -.2 .2 1000 1000 1000  image.bmp ; photoqt image.bmp &

run9: 
	./julia_acc_s 0 -0.4 0.6  -2 2 -2 2  1000 1000 1000  image.bmp ; photoqt image.bmp &

runx:
	./julia_acc_s 0 -0.4 0.6 -0.181862 -0.181772 1.019085 1.019175 4000 4000 4000 image.bmp;


clean:
	@rm -rf $(OBJS) julia_acc_s *~ *.bak *.bmp
