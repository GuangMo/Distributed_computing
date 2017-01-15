
all:
	make -C julia-omp julia ;
	make -C julia-acc-d julia ;
	make -C julia-acc-s julia ;
	cp julia-omp/julia_omp julia_omp;
	cp julia-acc-d/julia_acc_d julia_acc_d;
	cp julia-acc-s/julia_acc_s julia_acc_s;


clean:
	@rm -rf $(OBJS) julia_omp julia_acc_s julia_acc_d *~ *.bak *.bmp *.zip
	make -C julia-omp clean;
	make -C julia-acc-s clean;
	make -C julia-acc-d clean;

