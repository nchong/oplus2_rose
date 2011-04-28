
static void res(int *i, int *j, double **per_node_start, double *per_node, double *per_edge) {

  if (*j == 8) {
    printf("(2) per_node_start[0] = %p\n",  per_node_start[0]);
    printf("    per_node          = %p\n",  per_node);
    printf("    offset            = %ld\n", per_node - per_node_start[0]);
  }

  per_node[0] += 1;
  per_node[1] += 1;
  per_node[2] += 1;

  per_edge[0] += 1;
  per_edge[1] += 1;
  per_edge[2] += 1;

}
