#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <time.h>

int sub_nz;
int sub_ny;
int sub_nx;
int rank, size, PX, PY, PZ, NX, NY, NZ, NC;

int get(int nx,int ny,int nz){
    if(nx < 0){
        return -1;
    }
    if(ny < 0){
        return -1;
    }
    if(nz < 0){
        return -1;
    }
    int rank = nz*PY*PX + ny*PX * nx;
    if(rank >= 0 && rank < size){
        return rank;
    }
    return -1;
}

void transfer_function(int r1, int x,int y,int z, float *send_buf, float *recv_buf, int buf_size) {
    int r2 = get(x,y,z);
    if(r2 == -1){
        return;
    }
    MPI_Request reqs[2];
    MPI_Status stats[2];

    MPI_Isend(send_buf, buf_size, MPI_FLOAT, r2, r2, MPI_COMM_WORLD, &reqs[0]);
    MPI_Irecv(recv_buf, buf_size, MPI_FLOAT, r2, r1, MPI_COMM_WORLD, &reqs[1]);

    MPI_Waitall(2, reqs, stats);
}

int main(int argc, char **argv) {
    char input_file[256], output_file[256];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc != 10) {
        if (rank == 0) {
            printf("Usage: mpirun -np P ./executable <input_file> PX PY PZ NX NY NZ NC <output_file>\n");
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    sscanf(argv[1], "%s", input_file);
    PX = atoi(argv[2]);
    PY = atoi(argv[3]);
    PZ = atoi(argv[4]);
    NX = atoi(argv[5]);
    NY = atoi(argv[6]);
    NZ = atoi(argv[7]);
    NC = atoi(argv[8]);
    sscanf(argv[9], "%s", output_file);

    if (size != PX * PY * PZ) {
        if (rank == 0) {
            printf("Error: The number of processes must be equal to PX * PY * PZ.\n");
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int x_rank = rank % PX;
    int y_rank = (rank / PX) % PY;
    int z_rank = rank / (PX * PY);

    sub_nx = NX / PX;
    sub_ny = NY / PY;
    sub_nz = NZ / PZ;

    float *data = NULL;
    float *data1 = (float *)malloc(sub_nz * NY * NX * NC * sizeof(float));

    if (rank == 0) {
        data = (float *)malloc(NZ * NY * NX * NC * sizeof(float));

        FILE *file = fopen(input_file, "r");
        if (!file) {
            perror("Error opening file");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        for (int i = 0; i < NZ * NY * NX * NC; i++) {
            if (fscanf(file, "%f", &data[i]) != 1) {
                perror("Error reading file");
                fclose(file);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }
        fclose(file);
    }

    MPI_Comm z_comm;
    int cz = (x_rank == 0 && y_rank == 0) ? 1 : MPI_UNDEFINED;
    MPI_Comm_split(MPI_COMM_WORLD, cz, rank, &z_comm);
    if (cz == 1) {
        MPI_Scatter(data, sub_nz * NY * NX * NC, MPI_FLOAT, data1, sub_nz * NY * NX * NC, MPI_FLOAT, 0, z_comm);
        MPI_Comm_free(&z_comm);
    }
    if (rank == 0) free(data);

    float *data2 = (float *)malloc(sub_nz * sub_ny * NX * NC * sizeof(float));

    MPI_Comm n_comm;
    MPI_Comm_split(MPI_COMM_WORLD, z_rank, rank, &n_comm);

    MPI_Comm y_comm;
    int cy = (x_rank == 0) ? 1 : MPI_UNDEFINED;
    MPI_Comm_split(n_comm, cy, rank, &y_comm);

    if (cy == 1) {
        for (int i = 0; i < sub_nz; i++) {
            MPI_Scatter(&data1[i * NY * NX * NC], sub_ny * NX * NC, MPI_FLOAT, &data2[i * sub_ny * NX * NC], sub_ny * NX * NC, MPI_FLOAT, 0, y_comm);
        }
        MPI_Comm_free(&y_comm);
    }

    free(data1);

    float *data3 = (float *)malloc(sub_nz * sub_ny * sub_nx * NC * sizeof(float));

    MPI_Comm m_comm;
    MPI_Comm_split(n_comm, y_rank, rank, &m_comm);

    for (int k = 0; k < sub_nz; k++) {
        for (int j = 0; j < sub_ny; j++) {
            MPI_Scatter(&data2[(k * sub_ny + j) * NX * NC], sub_nx * NC, MPI_FLOAT, &data3[(k * sub_ny + j) * sub_nx * NC], sub_nx * NC, MPI_FLOAT, 0, m_comm);
        }
    }

    free(data2);
    MPI_Comm_free(&m_comm);
    MPI_Comm_free(&n_comm); 

    float data4[sub_nz][sub_ny][sub_nx][NC];
    float mn_val[sub_nz][sub_ny][sub_nx][NC];
    float mx_val[sub_nz][sub_ny][sub_nx][NC];
    float min_val = INT32_MAX, max_val = INT32_MIN;
    int min_count = 0;
    int max_count = 0;

    int index = 0;
    for(int i = 0;i<sub_nz;i++){
        for(int j = 0;j<sub_ny;j++){
            for(int k=0;k<sub_nx;k++){
                for(int x = 0;x < NC;x++){
                    data4[i][j][k][x] = data3[index];
                    mn_val[i][j][k][x] = INT32_MAX;
                    mx_val[i][j][k][x] = INT32_MIN;
                    index++;
                }
            }
        }
    }

    //data for faces
    float right[sub_nz][sub_ny][NC];
    float left[sub_nz][sub_ny][NC];
    float bottom[sub_ny][sub_nx][NC];
    float top[sub_ny][sub_nx][NC];
    float front[sub_nz][sub_nx][NC];
    float back[sub_nz][sub_nx][NC];
    for(int j = 0;j<sub_ny;j++){
        for(int i = 0;i<sub_nx;i++){
            for(int x = 0;x < NC ;x++){
                top[j][i][x] = data4[0][j][i][x];
                bottom[j][i][x] = data4[sub_nz-1][j][i][x];
            }
        }
    }
    for(int j = 0;j<sub_nz;j++){
        for(int i = 0;i<sub_nx;i++){
            for(int x = 0;x < NC ;x++){
                front[j][i][x] = data4[j][0][i][x];
                back[j][i][x] = data4[j][sub_ny-1][i][x];
            }
        }
    }
    for(int j = 0;j<sub_nz;j++){
        for(int i = 0;i<sub_ny;i++){
            for(int x = 0;x < NC ;x++){
                right[j][i][x] = data4[j][i][sub_nx-1][x];
                left[j][i][x] = data4[j][i][0][x];
            }
        }
    }

    //faces data gathering
    float recv_right[sub_nz][sub_ny][NC];
    float recv_left[sub_nz][sub_ny][NC];
    float recv_bottom[sub_ny][sub_nx][NC];
    float recv_top[sub_ny][sub_nx][NC];
    float recv_front[sub_nz][sub_nx][NC];
    float recv_back[sub_nz][sub_nx][NC];
    transfer_function(rank,x_rank,y_rank,z_rank -1,top,recv_top,sub_ny*sub_nx*NC);
    transfer_function(rank,x_rank,y_rank,z_rank+1,bottom,recv_bottom,sub_ny*sub_nx*NC);
    transfer_function(rank,x_rank-1,y_rank,z_rank,left,recv_left,sub_ny*sub_nz*NC);
    transfer_function(rank,x_rank+1,y_rank,z_rank,right,recv_right,sub_ny*sub_nz*NC);
    transfer_function(rank,x_rank,y_rank-1,z_rank,front,recv_front,sub_nx*sub_nz*NC);
    transfer_function(rank,x_rank,y_rank+1,z_rank,back,recv_back,sub_nx*sub_nz*NC);

    //data for edges
    float tf[sub_nx][NC],tb[sub_nx][NC],tr[sub_ny][NC],tl[sub_ny][NC];
    float bf[sub_nx][NC],bb[sub_nx][NC],br[sub_ny][NC],bl[sub_ny][NC];
    float fr[sub_nz][NC],fl[sub_nz][NC],bar[sub_nz][NC],bal[sub_nz][NC];
    for(int i = 0;i<sub_nx;i++){
        for(int x = 0;x < NC;x++){
            tf[i][x] = data4[0][0][i][x];
            tb[i][x] = data4[0][sub_ny-1][i][x];
            bf[i][x] = data4[sub_nz-1][0][i][x];
            bb[i][x] = data4[sub_nz-1][sub_ny-1][i][x];
        }
    }
    for(int i = 0;i<sub_ny;i++){
        for(int x = 0;x < NC;x++){
            tl[i][x] = data4[0][i][0][x];
            tr[i][x] = data4[0][i][sub_nx-1][x];
            bl[i][x] = data4[sub_nz-1][i][0][x];
            br[i][x] = data4[sub_nz-1][i][sub_nx-1][x];
        }
    }
    for(int i = 0;i<sub_nz;i++){
        for(int x = 0;x < NC;x++){
            fr[i][x] = data4[i][0][sub_nx-1][x];
            fl[i][x] = data4[i][0][0][x];
            bal[i][x] = data4[i][sub_ny-1][0][x];
            bar[i][x] = data4[i][sub_ny-1][sub_nx-1][x];
        }
    }

    //data gathering for edges ..
    float recv_tf[sub_nx][NC],recv_tb[sub_nx][NC],recv_tr[sub_ny][NC],recv_tl[sub_ny][NC];
    float recv_bf[sub_nx][NC],recv_bb[sub_nx][NC],recv_br[sub_ny][NC],recv_bl[sub_ny][NC];
    float recv_fr[sub_nz][NC],recv_fl[sub_nz][NC],recv_bar[sub_nz][NC],recv_bal[sub_nz][NC];

    transfer_function(rank,x_rank-1,y_rank,z_rank-1,tl,recv_tl,sub_nx*NC);
    transfer_function(rank,x_rank+1,y_rank,z_rank-1,tr,recv_tr,sub_nx*NC);
    transfer_function(rank,x_rank,y_rank+1,z_rank-1,tf,recv_tf,sub_ny*NC);
    transfer_function(rank,x_rank,y_rank-1,z_rank-1,tb,recv_tb,sub_ny*NC);
    transfer_function(rank,x_rank-1,y_rank,z_rank+1,bl,recv_bl,sub_nx*NC);
    transfer_function(rank,x_rank+1,y_rank,z_rank+1,br,recv_br,sub_nx*NC);
    transfer_function(rank,x_rank,y_rank+1,z_rank+1,bf,recv_bf,sub_ny*NC);
    transfer_function(rank,x_rank,y_rank-1,z_rank+1,bb,recv_bb,sub_ny*NC);
    transfer_function(rank,x_rank-1,y_rank-1,z_rank,fl,recv_fl,sub_nz*NC);
    transfer_function(rank,x_rank+1,y_rank-1,z_rank,fr,recv_fr,sub_ny*NC);
    transfer_function(rank,x_rank-1,y_rank+1,z_rank,bal,recv_bal,sub_ny*NC);
    transfer_function(rank,x_rank-1,y_rank+1,z_rank,bar,recv_bar,sub_ny*NC);

    //corners
    
    float tfr[NC],tfl[NC],tbr[NC],tbl[NC],bfr[NC],bfl[NC],bbr[NC],bbl[NC];
    for(int x =0 ;x<NC;x++){
        tfl[x] = data4[0][0][0][x];
        tfr[x] = data4[0][0][sub_nx-1][x];
        tbl[x] = data4[0][sub_ny-1][0][x];
        tbr[x] = data4[0][sub_ny-1][sub_nx-1][x];
        bfl[x] = data4[sub_nz-1][0][0][x];
        bfr[x] = data4[sub_nz-1][0][sub_nx-1][x];
        bbl[x] = data4[sub_nz-1][sub_ny-1][0][x];
        bbr[x] = data4[sub_nz-1][sub_ny-1][sub_nx-1][x];
    }

    
    float recv_tfr[NC],recv_tfl[NC],recv_tbr[NC],recv_tbl[NC],recv_bfr[NC],recv_bfl[NC],recv_bbr[NC],recv_bbl[NC];

    transfer_function(rank,x_rank+1,y_rank-1,z_rank-1,tfr,recv_tfr,NC);
    transfer_function(rank,x_rank-1,y_rank-1,z_rank-1,tfl,recv_tfl,NC);
    transfer_function(rank,x_rank+1,y_rank+1,z_rank-1,tbr,recv_tbr,NC);
    transfer_function(rank,x_rank-1,y_rank+1,z_rank-1,tbl,recv_tbl,NC);
    transfer_function(rank,x_rank+1,y_rank-1,z_rank+1,bfr,recv_bfr,NC);
    transfer_function(rank,x_rank-1,y_rank-1,z_rank+1,bfl,recv_bfl,NC);
    transfer_function(rank,x_rank+1,y_rank+1,z_rank+1,bbr,recv_bbr,NC);
    transfer_function(rank,x_rank-1,y_rank+1,z_rank+1,bbl,recv_bbl,NC);

    

    

    

    free(data3);



    MPI_Finalize();
    return 0;
}
