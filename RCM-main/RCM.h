typedef struct V1
{
    int ii;
    struct V1* next;

} VERTEX;


typedef struct
{
    int n;

    int* R;
    VERTEX* Q;
    VERTEX* Qfinal;    
    
    int* deg;
    
    VERTEX** nei;

    bool* Rf;
    bool* Qf;    
    
} RCM;

RCM* rcmInit(MESH* mesh);

void rcmCalcNei(RCM* rcm, MESH* mesh);

int rcmS1(RCM* rcm, MESH* mesh);

void rcmSortNei(RCM* rcm, int index);

void rcmPushQ(RCM* rcm, VERTEX* C);
