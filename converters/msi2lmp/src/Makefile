
TARGET	 = msi2lmp.exe

SRCS     = msi2lmp.c \
           ReadCarFile.c \
           ReadMdfFile.c \
	   MakeLists.c \
           ReadFrcFile.c \
           InitializeItems.c \
	   SearchAndFill.c \
           GetParameters.c \
	   CheckLists.c \
           WriteDataFile.c

OBJS     = $(SRCS:.c=.o)

HEADERS  = Msi2LMP2.h Forcefield.h

CC       = gcc
CFLAGS   = -O2
FRCFILE  = cvff.frc
FRCFILE2 = cff91.frc
README   = README
MKFILE   = Makefile

$(TARGET) : $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) -lm

$(OBJS) : $(HEADERS)

.c.o:
	$(CC) $(CFLAGS) -c $<
clean:
	rm -f $(OBJS) $(TARGET)


