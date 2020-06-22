
Il makefile prevede le seguenti opzioni:
make: compila ed esegue il codice con input.dat
make generalEq: equilibra in maniera automatica con qualche run da 1000 step
make Argon"fase"Eq copia il file di input giusto che permette di effettuare l'equilibrazione del sistema e di conseguenza, grazie alla prima riga del file di input, di salvare i file direttamente nelle cartelle della fase dell'argon giusta.
make clean: rimuove le precedenti configurazioni, necessario prima di fare l'equilibrazione per ripartire dalla configurazione iniziale
L'algoritmo prevede la possibilità di riscalare le velocità e di riprendere da una configurazione di una run precedente: ciò è fatto tramite due interi 1/0 nel file di input (rescaling e restart)
