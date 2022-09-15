double **mat_adj(double **points, int dim, int len){
    int i, j;
    double temp;
    double **a = alloc_mat(len, len);
    if (a == NULL)
        return NULL;
    i = 0;
    while (i < len)
    {
        j = 0;
        while (j < len)
        {
           if(i == j)
               a[i][j] = 0;
            else {
                temp = (exp((euc_norm_calc(points[i], points[j], dim)) / (-2)));
                a[i][j] =temp; 
            }
           j++;
        }   
        i++;
    }
    return a;
}


double euc_norm_calc(double *arr1, double *arr2, int dim){
    int j;
    double sum, temp;
    sum = 0;
    j = 0;
    while(j < dim){
        sum = pow(arr1[j] - arr2[j], 2) + sum;
        j++;
    }
    temp = sqrt(sum);
    sum = temp;
    return sum;
}
