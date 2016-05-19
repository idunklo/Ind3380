 DEBUG ?= OFF
 ifeq ($(DEBUG), ON)
		CFLAGS = -g -std=c99 -openmp -DDEBUG
else
		CFLAGS = -std=c99 -openmp -O3
endif

CC = mpicc
LDFLAGS = -lm

EXEC = oblig2
OBJS = oblig2.o

$(EXEC) : $(OBJS)
		$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

obj-clean:
		$(RM) *.o

exec-clean:
		$(RM) $(EXEC)

clean:
		$(MAKE) obj-clean
			$(MAKE) exec-clean
