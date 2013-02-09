#include <sys/types.h>
#include <signal.h>
#include <unistd.h>
 
inline static void core_dump(int sigid)
{
    abort();
}
