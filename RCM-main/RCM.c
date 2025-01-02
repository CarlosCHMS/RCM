#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include <stdbool.h>
#include<math.h>
#include"utils.h"
#include"mesh.h"
#include"RCM.h"

RCM* rcmInit(MESH* mesh)
{
    RCM* rcm = malloc(sizeof(RCM));

    rcm->R = malloc(mesh->Nelem*sizeof(int));
    rcm->deg = malloc(mesh->Nelem*sizeof(int));    
    rcm->Q = NULL;
    rcm->n = 0;

    rcm->Rf = malloc(mesh->Nelem*sizeof(bool));
    rcm->Qf = malloc(mesh->Nelem*sizeof(bool));

    rcm->nei = malloc(mesh->Nelem*sizeof(VERTEX*));
    
    for(int ii=0; ii<mesh->Nelem; ii++)
    {
        rcm->Rf[ii] = false;
        rcm->Qf[ii] = false;        
        rcm->nei[ii] = NULL;
        rcm->deg[ii] = 0;
    }
    
    rcmCalcNei(rcm, mesh);

    return rcm;
}


void rcmCalcNei(RCM* rcm, MESH* mesh)
{
    VERTEX* v;
    int e0, e1;    

    for(int ii=0; ii<mesh->Ncon; ii++)
    {
        e0 = mesh->con[ii][0];
        e1 = mesh->con[ii][1];
        
        rcm->deg[e0] += 1;
        rcm->deg[e1] += 1;        

        v = malloc(sizeof(VERTEX));
        v->ii = e1;
        v->next = rcm->nei[e0];
        rcm->nei[e0] = v;

        v = malloc(sizeof(VERTEX));
        v->ii = e0;
        v->next = rcm->nei[e1];
        rcm->nei[e1] = v;        
    }
    
    for(int ii=0; ii<mesh->Nelem; ii++)
    {
        rcmSortNei(rcm, ii);
    }
    
    /*
    for(int ii=0; ii<mesh->Nelem; ii++)
    {
        v = rcm->nei[ii];
        while(v != NULL)
        {
            printf("%i,", rcm->deg[v->ii]);
            v = v->next;
        }
        printf("\n");
    }
    */
}


int rcmS1(RCM* rcm, MESH* mesh)
{
    int deg = rcm->deg[0];
    int index = 0;
    
    for(int ii=1; ii<mesh->Nelem; ii++)
    {
        if(rcm->deg[ii] < deg)
        {
            deg = rcm->deg[ii];
            index = ii;
        }
    }
    
    return index;
}


void rcmSortNei(RCM* rcm, int index)
{
    int iaux;
    bool flag = true;
    VERTEX* v0;
    VERTEX* v1;

    while(flag)
    {
        flag = false;
        v0 = rcm->nei[index];
        v1 = v0->next;
        while(v1 != NULL)
        {            
            if(rcm->deg[v0->ii] > rcm->deg[v1->ii])
            {
                iaux = v0->ii;
                v0->ii = v1->ii;
                v1->ii = iaux;
                flag = true;
            }
          
            v0 = v1;
            v1 = v1->next;
        }
    }
}

void rcmPushQ(RCM* rcm, VERTEX* C)
{
    if(rcm->Q == NULL)
    {
        C->next = NULL;
        rcm->Q = C;
        rcm->Qfinal = C;
    }
    else
    {
        C->next = NULL;
        rcm->Qfinal->next = C;
        rcm->Qfinal = C;
    }
}
