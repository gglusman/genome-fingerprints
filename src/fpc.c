/*
 * fpc.c
 *
 *  Created on: Jul 5, 2017
 *      Author: mrobinso
 */

#define _POSIX_SOURCE

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#ifndef TRUE
#  define TRUE (0 == 0)
#endif /* TRUE */

#ifndef FALSE
#  define FALSE (! TRUE)
#endif /* FALSE */

#define ONE ((size_t) 1)

float minValue = -2.0;

typedef struct _FingerprintSet {
    int  fingerprints;
    int  length;
    int *valuep;
} FingerprintSet;


/*--------------------------------------------------------------------------------
// readFingerprintSet("queryFile")
*/
static void
readFingerprintSet(char *filename, FingerprintSet *fpp)
{
    FILE  *fp   = fopen(filename, "r");
    size_t size = (size_t) 0;
    int length;

    if ((FILE *) NULL != fp) {
        if (ONE == fread((void *) &length, sizeof(int), ONE, fp)) {
            fpp->length = length;
            fseek(fp, 0, SEEK_END);
            size = ftell(fp);
            fpp->fingerprints = (int) (size / sizeof(int)) - 1;
            if (0 == (fpp->fingerprints % fpp->length)) {
                fpp->fingerprints /= fpp->length;
                fpp->valuep = (int *)
                    calloc(fpp->fingerprints, fpp->length*sizeof(int));
                fseek(fp, (unsigned long int) sizeof(int), SEEK_SET);
                size_t k = fread((void *) fpp->valuep,
                                 fpp->length*sizeof(int),
                                 (size_t) fpp->fingerprints, fp);
                if ((size_t) fpp->fingerprints != k) {
                    fprintf(stderr, "Couldn't read %lu fingerprints from file '%s'\n",
                            (unsigned long) fpp->fingerprints, filename);
                    fflush(stderr);
                    exit(0);
                }
                fclose(fp);
            } else {
                fprintf(stderr,
                    "ERROR: File '%s' doesn't contain fingerprints of size %ld\n",
                    filename, fpp->length);
                fflush(stderr);
                exit(0);
            }
        }
    } else {
        fprintf(stderr,"ERROR: Couldn't read from fingerprint file '%s'\n", filename);
        fflush(stderr);
        exit(0);
    }
}

/*--------------------------------------------------------------------------------
// compare(queryp, targetp)
//--------------------------------------------------------------------------------
*/

static void
selfCompare(FingerprintSet *queryp)
{
    int   *qEndp = queryp->valuep + (queryp->fingerprints  * queryp->length);
    float *qBufp, *tBufp, *tVarp;
    size_t alloc = 2*queryp->length + queryp->fingerprints;
    size_t bufSize = sizeof(float)*queryp->length;

    /* Allocate working space */
    qBufp = (float *) calloc(sizeof(float), alloc);
    if ((float *) NULL == qBufp) {
        fprintf(stderr, "ERROR: Couldn't allocate %lu floats\n",
                (unsigned long int) alloc);
        fflush(stderr);
        exit(0);
    }
    tVarp = (tBufp = qBufp + queryp->length) + queryp->length;

    {
        int   query = 0;
        float n     = (float) queryp->length;
        float qVar, tVar;
        register int  *qp;

        /* Center all queries */
        for (qp = queryp->valuep; qp < qEndp; qp += queryp->length) {
            register float  mean = 0.0;
            register int   *xp = qp;
            register float *qbp = qBufp;
            register float *tbp;

            while (qbp < tBufp) {
                mean += ((*qbp++) = (float) *xp++);
            }
            mean /= n;

            /* Center the query, compute its variance */
            qVar = 0.0;
            for (qbp = qBufp; qbp < tBufp; ++qbp) {
                (*qbp) -= mean;
                qVar   += (*qbp)*(*qbp);
            }
            /* Save the centered query */
            (void) memcpy((void *) qp, (void *) qBufp, bufSize);
            tVarp[query++] = (float) sqrt((double) qVar);
        }

        /* Perform the comparisons */
        query = 0;
        for (qp = queryp->valuep; qp < qEndp; qp += queryp->length) {
            register float *qbp;
            register float *tbp;
            register int   *tp;

            /* Recall the query */
            qVar = tVarp[query++];
            (void) memcpy((void *) qBufp, (void *) qp, bufSize);

            int target = query;
            for (tp  = qp + queryp->length; tp < qEndp; tp += queryp->length) {
                register float r = 0.0;
                /* Recall the target */
                tVar = tVarp[target];
                (void) memcpy((void *) tBufp, (void *) tp, bufSize);

                qbp = qBufp;
                tbp = tBufp;
                while (tbp < tVarp) {
                    r += (*qbp++) * (*tbp++);
                }
                r /= (qVar*tVar);
                ++ target;
                if (minValue <= r) {
                    printf("%d\t%d\t%7f\n", query, target, r);
                }
            }
        }
    }
    free((void *) qBufp);
}


static void
compare1(FingerprintSet *queryp,
         FingerprintSet *targetp,
         int             flipped)
{
    int *qEndp = queryp->valuep  + (queryp->fingerprints  * queryp->length);
    int *tEndp = targetp->valuep + (targetp->fingerprints * targetp->length);
    float    *qBufp, *tBufp, *tVarp;
    if (queryp->length == targetp->length) {

        /* Allocate working space */
        size_t bufSize = sizeof(float)*targetp->length;
        size_t alloc   = 2*queryp->length;
        if (1 < queryp->fingerprints) alloc += targetp->fingerprints;
        qBufp = (float *) calloc(sizeof(float), alloc);
        if ((float *) NULL == qBufp) {
            fprintf(stderr, "ERROR: Couldn't allocate %lu floats\n",
                    (unsigned long int) alloc);
            fflush(stderr);
            exit(0);
        }
        tVarp = (tBufp = qBufp + queryp->length) + queryp->length;

        /* First query */
        register float  mean = 0.0;
        register int   *qp   = queryp->valuep;
        register int   *xp   = qp;
        register float *qbp  = qBufp;
        float n = (float) queryp->length;

        /* Compute the query mean */
        while (qbp < tBufp) {
            mean += ((*qbp++) = (float) *xp++);
        }
        mean /= n;

        /* Center the query, compute its variance */
        float qVar = 0.0;
        for (qbp = qBufp; qbp < tBufp; ++qbp) {
            (*qbp) -= mean;
            qVar   += (*qbp)*(*qbp);
        }
        qVar = (float) sqrt((double) qVar);

        /* Process each target */
        int  query  = 1;
        int  target = 0;
        int *tp;
        for (tp  = targetp->valuep;
             tp < tEndp; tp += targetp->length) {
            register float *tbp = tBufp;

            /* Compute the target mean */
            mean = 0.0;
            xp   = tp;
            while (tbp < tVarp) {
                mean += ((*tbp++) = (float) *xp++);
            }
            mean /= n;

            /* Center the target; compute its variance; compute Pearson's r */
            float tVar = 0.0;
            register float r = 0.0;
            qbp = qBufp;
            for (tbp = tBufp; tbp < tVarp; ++tbp) {
                (*tbp) -= mean;
                tVar   += (*tbp)*(*tbp);
                r      += (*tbp)*(*qbp++);
            }
            tVar = (float) sqrt((double)tVar);

            if (1 < queryp->fingerprints) { /* If there is only one query, skip this */
                tVarp[target] = tVar;
                (void) memcpy((void *) tp, (void *) tBufp, bufSize);
            }

            r /= qVar*tVar;
            ++ target;

            if (minValue <= r) {
                if (flipped) {
                    printf("%d\t%d\t%7f\n", target, query, r);
                } else {
                    printf("%d\t%d\t%7f\n", query, target, r);
                }
            }
        }

        /* Subsequent queries */
        for (qp += queryp->length; qp < qEndp; qp += queryp->length) {
            ++ query;

            /* Compute the query mean */
            mean = 0.0;
            xp   = qp;
            qbp  = qBufp;
            while (qbp < tBufp) {
                mean += ((*qbp++) = (float) *xp++);
            }
            mean /= n;

            /* Center the query, compute its variance */
            for (qbp = qBufp; qbp < tBufp; ++qbp) {
                (*qbp) -= mean;
                qVar   += (*qbp)*(*qbp);
            }
            qVar = (float) sqrt((double) qVar);

            target = 0;
            for (tp  = targetp->valuep;
                 tp < tEndp; tp += targetp->length) {
                register float  r = 0.0;
                register float *tbp = tBufp;

                /* Recall the target */
                float tVar = tVarp[target];
                (void) memcpy((void *) tBufp, (void *) tp, bufSize);

                qbp = qBufp;
                while (tbp < tVarp) {
                    r += (*qbp++)*(*tbp++);
                }
                r /= (qVar*tVar);
                ++ target;

                if (flipped) {
                    printf("%d\t%d\t%7f\n", target, query, r);
                } else {
                    printf("%d\t%d\t%7f\n", query, target, r);
                }
            }
        }
        free((void *) qBufp);
    } else {
        fprintf(stderr, "ERROR: Fingerprints are of different lengths (%ld vs %ld)\n",
                queryp->length, targetp->length);
        fflush(stderr);
        exit(0);
    }
}


static void
compare(FingerprintSet *queryp,
        FingerprintSet *targetp)
{
    if (targetp->fingerprints < queryp->fingerprints) {
        compare1(targetp, queryp, (int) (0 == 0)); /* TRUE */
    } else {
        compare1(queryp, targetp, (int) (0 != 0)); /* FALSE */
    }
}

/*--------------------------------------------------------------------------------
//
*/
static char *USAGE = "Usage: %s [ -t <threshold> ] <queryFile> [ <dbFile1> ... ]\n";

int
main(int   argc,
     char *argv[])
{
    FingerprintSet query;
    extern int   optind;
    extern char *optarg;
    int opt;
    int usage = FALSE;

    while (-1 != (opt = getopt(argc,argv,"t:"))) {
        switch (opt) {
        case 't':
            if (1 != sscanf(optarg,"%f",&minValue)) {
                usage = TRUE;
            }
            break;
        default:
            usage = TRUE;
            break;
        }
    }
    if (argc == optind) { usage = TRUE; } 
    if (! usage) {
        readFingerprintSet(argv[optind++],&query);
        if (optind == argc) {
            selfCompare(&query);
        } else {
            FingerprintSet db;
            while (optind < argc) {
                readFingerprintSet(argv[optind++],&db);
                compare(&query, &db);
                free((void *) db.valuep);
            }
        }
        free((void *) query.valuep);
    }
    if (usage) {
        fprintf(stderr, USAGE, argv[0]);
        fflush(stderr);
    }
}

/* end of fpcompare.c */
