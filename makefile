all:
	make -C CTransR
	make -C PTransE
	make -C TransE
	make -C TransH
	make -C TransR
	make -C cluster

