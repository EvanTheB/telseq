telseq_debug: telseq.c
	gcc -Wno-unknown-pragmas -Wno-unused-parameter -I/home/evaben/download/htslib -Wall -Wextra -g -O2 -std=c99 -o telseq.o -c telseq.c
	gcc -Wall -Wextra -g -O2 -o telseq_debug telseq.o /home/evaben/download/htslib/libhts.a -lz -lcurl /home/evaben/download/libdeflate/libdeflate.a -lz -lcrypto -lpthread -llzma -lbz2

clean:
	rm telseq_debug
