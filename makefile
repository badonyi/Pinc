CC ?= gcc
CFLAGS ?= -O3 -Wall -Wextra
LDLIBS = -lm

ifeq ($(OS),Windows_NT)
	TARGET = Pinc.exe
	CLEAN_CMD = del /Q /F $(TARGET)
else
	TARGET = Pinc
	CLEAN_CMD = rm -f $(TARGET)
endif

$(TARGET): src/pinc.c
	$(CC) $(CFLAGS) src/pinc.c -o $(TARGET) $(LDLIBS)

clean:
	$(CLEAN_CMD)

.PHONY: clean
