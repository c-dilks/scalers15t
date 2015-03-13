
all: scaler2_reader_bit.exe

scaler2_reader_bit.exe :	scaler2_reader_bit.c
	gcc scaler2_reader_bit.c -o scaler2_reader_bit.exe

clean:	
	/bin/rm -f *_bit.o *~ scaler2_reader_bit.exe
