CC     = gcc
CFLAGS = -march=native -O2 -Wall

SRC = convert_binary_output_file_to_ascii.c
EXE = convert_binary_output_file_to_ascii

all: $(EXE)

$(EXE): $(SRC)
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f $(EXE)
