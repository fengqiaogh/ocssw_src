#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stdint.h>

#define isdigit(c) (c >= '0' && c <= '9')
#define isalpha(c) ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z'))
#define isspace(c)   (c == ' ' || c == '\t' || c == '\n')

int zs_strparts(char str[], char delim[], int maxparts, char *strvec[]);
int zs_strip(char str[]);
char *zs_substr(char *longstr, char *smlstr);

/**
 * This function parses a path name.  It removes extra /'s,
 * translates any environment variables,
 * translates ~users, and also handles "." and ".."
 *
 * @param inpath file name to be translated
 * @param outpath pointer to enough space for the translated file path
 */
void parse_file_name(const char* inpath, char* outpath) {
    char str1[FILENAME_MAX], str2[FILENAME_MAX], envvar[128];
    char *dollar, *ptr, *trans, *parts[128];
    int32_t numparts, i, index1, index2;

    strcpy(outpath, "");
    /*
         First translate all environment variable3s in the input string
     */
    ptr = (char*) &inpath[0];
    index1 = 0;
    strcpy(str1, "");

    while ((dollar = strchr(ptr, '$'))) {
        strncat(str1, ptr, (dollar - ptr));
        index1 += (dollar - ptr);
        str1[index1] = '\0';

        index2 = 0;
        ptr = dollar + 1;
        while (isalpha(*ptr) || isdigit(*ptr) || *ptr == '_') {
            envvar[index2++] = *ptr;
            ptr++;
        }
        envvar[index2] = '\0';
        /*
           Translate the variable
         */
        if ((trans = getenv(envvar)) == NULL) /* Environment variable Not Found */ {
            strcat(str1, "$");
            strcat(str1, envvar);
            index1 += (strlen(envvar) + 1);
        } else {
            strcat(str1, trans);
            index1 += strlen(trans);
        }
    }
    strcat(str1, ptr);
    /*
       Now all environment variables are translated.
       Break up the string into parts separated by /
     */
    numparts = zs_strparts(str1, "/", 128, parts);

    for (i = 0; i < numparts; i++)
        zs_strip(parts[i]);
    index1 = 0;
    strcpy(str2, "");

    if (strlen(parts[0]) == 0) /* either empty or first char is '/' */ {
        for (i = 0; i < numparts; i++) {
            if ((strlen(parts[i]) == 0) || (strcmp(parts[i], ".") == 0))
                index1++;
            else
                break;
        }
        if ((i == numparts) && (numparts != 1))
            strcpy(str2, "/");
        else if ((i < numparts) && (strcmp(parts[i], "..") == 0))
            strcpy(str2, "/");
    } else if (*parts[0] == '~') /* get a home directory */ {
        printf("\nInterpretation of ~ is no longer supported.\n");
        exit(1);
    }
    else if ((strcmp(parts[0], ".") == 0) || (strcmp(parts[0], "..") == 0)) {
        getcwd(str2, 257);
        if (strcmp(parts[0], ".") == 0)
            index1 = 1;
    }

    for (i = index1; i < numparts; i++) {
        if (strcmp(parts[i], "..") == 0) {
            if ((ptr = strrchr(str2, '/')) == NULL)
                strcpy(str2, "");
            else {
                if (ptr != &str2[0])
                    *ptr = '\0';
                else
                    strcpy(str2, "/");
            }
        } else if (strlen(parts[i]) && (strcmp(parts[i], ".") != 0)) {
            if ((i != 0) && (strcmp(str2, "/") != 0))
                strcat(str2, "/");
            strcat(str2, parts[i]);
        }
    }

    strcpy(outpath, str2);
    for (i = 0; i < numparts; i++)
        free(parts[i]);
    return;
}

int zs_strparts(str, delim, maxparts, strvec)
char str[];
char delim[];
int maxparts;
char *strvec[];
{
    int i;
    int len;
    int numparts = 0;
    char *sptr1;
    char *sptr2;

    sptr1 = str;
    for (i = 0; i < maxparts; i++) {
        if (*sptr1 == '\0' || sptr1 == zs_substr(sptr1, delim))
            sptr2 = sptr1;
        else
            sptr2 = zs_substr(sptr1 + 1, delim);

        if (sptr2 == NULL)
            len = strlen(sptr1);
        else
            len = sptr2 - sptr1;

        strvec[i] = (char *) malloc(len + 1);
        strncpy(strvec[i], sptr1, len);
        *(strvec[i] + len) = '\0';

        numparts++;

        if (sptr2 == NULL || *sptr2 == '\0')
            break;
        sptr1 = sptr2 + 1;
    }

    return (numparts);
}

int zs_strip(str)
char str[];
{
    int32_t i, k, l;

    i = 0;
    l = strlen(str);
    if (isspace(str[i])) {
        while (isspace(str[i])) i++;
        //   strcpy (&str[0], &str[i]);
        // Change to memcpy()  JMG  02/13/12
        for (k = i; k <= l; k++)
            memcpy(&str[k - i], &str[k], 1);
    }

    i = strlen(str);
    while (--i >= 0)
        if (!isspace(str[i]))
            break;

    str[i + 1] = '\0';

    return (0);
}

char *zs_substr(longstr, smlstr)
char *longstr, *smlstr;
{
    char *pos;
    int llen,
            slen;

    slen = strlen(smlstr);
    llen = strlen(longstr);
    for (pos = longstr; *pos != '\0'; pos++) {
        if (slen > llen)
            return (NULL);
        if (!strncmp(pos, smlstr, slen))
            return (pos);
        llen--;
    }
    return (NULL);
}

