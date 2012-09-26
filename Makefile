# make is unsupported.  Use setup.py
CFLAGS=-pthread -g -DNDEBUG
LIBS=-lbam -lm

SOURCES=$(wildcard src/**/*.c src/*.c)
OBJECTS=$(patsubst %.c,%.o,$(SOURCES))

TARGET=build/tactmod.a
SO_TARGET=$(patsubst %.a,%.so,$(TARGET))

all: $(TARGET) $(SO_TARGET)

$(TARGET): CFLAGS += -fPIC
$(TARGET): build $(OBJECTS)
	ar rscu $@ $(OBJECTS)

$(SO_TARGET): $(TARGET) $(OBJECTS)
	$(CC) -shared -o $@ $(OBJECTS) 

build:
	mkdir -p build
clean:
	rm -f build/*
