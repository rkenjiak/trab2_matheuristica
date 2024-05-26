TRACE=NDEBUG
CFLAGS=-D NO_CONFIG_HEADER -D SCIP_VERSION_MAJOR -D$(TRACE) -g -std=c11


bin/mochila: bin/cmain.o bin/probdata_mochila.o bin/problem.o bin/heur_aleatoria.o bin/heur_rf.o bin/heur_lns.o bin/utils.o
	gcc $(CFLAGS) $(LDFLAGS) -o bin/mochila bin/cmain.o bin/probdata_mochila.o bin/problem.o  bin/heur_aleatoria.o  bin/heur_rf.o bin/heur_lns.o bin/utils.o -lscip -lm

bin/cmain.o: src/cmain.c
	gcc $(CFLAGS) -c -o bin/cmain.o src/cmain.c

bin/probdata_mochila.o: src/probdata_mochila.c src/probdata_mochila.h
	gcc $(CFLAGS) -c -o bin/probdata_mochila.o src/probdata_mochila.c

bin/problem.o: src/problem.c src/problem.h
	gcc $(CFLAGS) -c -o bin/problem.o src/problem.c

bin/heur_aleatoria.o: src/heur_aleatoria.c src/heur_aleatoria.h
	gcc $(CFLAGS) -c -o bin/heur_aleatoria.o src/heur_aleatoria.c

bin/heur_rf.o: src/heur_rf.c src/heur_rf.h
	gcc $(CFLAGS) -c -o bin/heur_rf.o src/heur_rf.c

bin/heur_lns.o: src/heur_lns.c src/heur_lns.h
	gcc $(CFLAGS) -c -o bin/heur_lns.o src/heur_lns.c

bin/utils.o: src/utils.c src/utils.h
	gcc $(CFLAGS) -c -o bin/utils.o src/utils.c

.PHONY: clean

clean:
	rm -f bin/*.o bin/mochila

