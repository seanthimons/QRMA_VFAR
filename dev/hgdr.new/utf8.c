
#include <stdio.h>
#include <string.h>

int grab_utf8(unsigned char *str, unsigned len, unsigned *target);

char buff[1000];
int main()
{
/* unsigned char pipo[] = { 'a', 'b' , 0xef ,0xbe ,0x88, 'x' ,'y' ,'z' , 0 }; */
unsigned idx ,val, len;
int ret;

while (fgets(buff, sizeof buff, stdin)) {
	len = strlen(buff);
	for(idx=0; idx < len;   ) {
		ret = grab_utf8((unsigned char*) buff+idx, len-idx, &val);
		/* fprintf(stderr, "%2x: %d: %x\n", idx, ret, val); */
		if (ret==1) {
			fprintf(stdout,"%c", val);
			}
		else    {
			fprintf(stdout,"# %x #", val);
			}
		if (ret>0) idx += ret;
		else idx++;
		}
	}
return 0;
}

int grab_utf8(unsigned char *str, unsigned len, unsigned *target)
{
unsigned idx;
unsigned val = 0;
unsigned todo;

if (!len) return 0;

val = str[0];
if ((val & 0x80) == 0x00) { *target = val; return 1; }
else if ((val & 0xe0) == 0xc0) { val &= 0x1f; todo = 1; }
else if ((val & 0xf0) == 0xe0) { val &= 0x0f; todo = 2; }
else if ((val & 0xf8) == 0xf0) { val &= 0x07; todo = 3; }
else if ((val & 0xfc) == 0xf8) { val &= 0x03; todo = 4; }
else if ((val & 0xfe) == 0xfc) { val &= 0x01; todo = 5; }
else {*target = val; return 1; } /* Default (Not in the spec) */
/* fprintf(stderr, "[val=%x, todo=%d]", val, todo); */

len--;str++;
if (todo > len) { return -todo; }

for(len=todo;todo--;) {
	val <<= 6;
	val |= *str++ & 0x3f;
	}

*target = val;
return 1+ len;
}

