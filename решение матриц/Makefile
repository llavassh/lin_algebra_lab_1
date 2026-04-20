CXX = g++
CXXFLAGS = -std=c++17 -O2 -Iinclude
TARGET = lab1.exe
SRCDIR = src
BUILDDIR = build

SOURCES = $(wildcard $(SRCDIR)/*.cpp)
OBJECTS = $(patsubst $(SRCDIR)/%.cpp, $(BUILDDIR)/%.o, $(SOURCES))

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $(TARGET)

$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp
	if not exist $(BUILDDIR) mkdir $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	if exist $(BUILDDIR) rmdir /s /q $(BUILDDIR)
	if exist $(TARGET) del $(TARGET)

run: $(TARGET)
	$(TARGET)

.PHONY: all clean run