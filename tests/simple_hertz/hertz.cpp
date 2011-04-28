#include "op_seq.h"
#include "res.h"
#include <fstream>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using namespace std;

// --------------------------------------------------------------------------
// UTILITIES
// --------------------------------------------------------------------------

//datastructure from serialized data (input and expected_output)
struct params {
  //constants
  double dt;
  double nktv2p;
  int ntype;
  double *yeff;
  double *geff;
  double *betaeff;
  double *coeffFrict;

  //node data
  int nnode;
  double *x;
  double *v;
  double *omega;
  double *radius;
  double *mass;
  int    *type;
  double *force;
  double *torque;

  //edge data
  int nedge;
  int *edge;
  double *shear;

  //partition data (OP2 only)
  int npartition;
  std::vector<int> partition_length;

  //expected results
  double *expected_force;
  double *expected_torque;
  double *expected_shear;
};

//unpickle array
template<class T>
inline void fill_array(std::ifstream &file, T *array, int num_elements) {
  if (file.eof()) {
    cout << "Error unexpected eof!" << endl;
    exit(-1);
  }
  for (int i=0; i<num_elements; i++) {
    file >> array[i];
  }
}

//unpickle file
struct params *parse_file(std::string fname) {
  ifstream file (fname.c_str(), ifstream::in);
  if (!file.is_open()) {
    cout << "Could not open [" << fname << "]" << endl;
    exit(-1);
  }
  if (file.bad()) {
    cout << "Error with file [" << fname << "]" << endl;
    exit(-1);
  }

  struct params *result = new params;
  int ntype;
  int nnode;
  int nedge;

  //constants
  file >> result->dt;
  file >> result->nktv2p;
  file >> ntype; result->ntype = ntype;
  result->yeff       = new double[ntype*ntype];
  result->geff       = new double[ntype*ntype];
  result->betaeff    = new double[ntype*ntype];
  result->coeffFrict = new double[ntype*ntype];
  fill_array(file, result->yeff,       (ntype*ntype));
  fill_array(file, result->geff,       (ntype*ntype));
  fill_array(file, result->betaeff,    (ntype*ntype));
  fill_array(file, result->coeffFrict, (ntype*ntype));

  //node data
  file >> nnode; result->nnode = nnode;
  result->x      = new double[nnode*3];
  result->v      = new double[nnode*3];
  result->omega  = new double[nnode*3];
  result->radius = new double[nnode  ];
  result->mass   = new double[nnode  ];
  result->type   = new int[nnode];
  result->force  = new double[nnode*3];
  result->torque = new double[nnode*3];
  fill_array(file, result->x,      nnode*3);
  fill_array(file, result->v,      nnode*3);
  fill_array(file, result->omega,  nnode*3);
  fill_array(file, result->radius, nnode);
  fill_array(file, result->mass,   nnode);
  fill_array(file, result->type,   nnode);
  fill_array(file, result->force,  nnode*3);
  fill_array(file, result->torque, nnode*3);

  //edge data
  file >> nedge; result->nedge = nedge;
  result->edge = new int[nedge*2];
  result->shear = new double[nedge*3];
  fill_array(file, result->edge,  nedge*2);
  fill_array(file, result->shear, nedge*3);

  //partition data (OP2 only)
  //by default no partition data means we just use the unpickled edge ordering
  result->npartition = 1;
  result->partition_length = vector<int>(1, nedge);

  //expected results
  result->expected_force  = new double[nnode*3];
  result->expected_torque = new double[nnode*3];
  result->expected_shear = new double[nedge*3];
  fill_array(file, result->expected_force,  nnode*3);
  fill_array(file, result->expected_torque, nnode*3);
  fill_array(file, result->expected_shear, nedge*3);

  return result;
}

// --------------------------------------------------------------------------
// MAIN
// --------------------------------------------------------------------------

int main(int argc, char **argv) {

  if (argc < 2) {
    printf("Usage: %s <stepfile>\n", argv[0]);
    return(1);
  }
  std::string step_filename(argv[1]);
  struct params *input = parse_file(step_filename);

  // HOST datastructs
  // unzip edge map so we know which node we're mutating in the kernel
  int *nodei = new int[input->nedge];
  int *nodej = new int[input->nedge];
  for (int e=0; e<input->nedge; e++) {
    nodei[e] = input->edge[(e*2)  ];
    nodej[e] = input->edge[(e*2)+1];
  }
  double *per_node = new double[input->nnode*3];
  for (int i=0; i<input->nnode*3; i++) {
    per_node[i] = 0.0;
  }
  double *per_edge = new double[input->nedge*3];
  for (int i=0; i<input->nedge*3; i++) {
    per_edge[i] = 0.0;
  }

  // OP2 init
  const char *fake_argv = {"foo"};
  op_init(1, (char **)fake_argv);

  // OP2 declarations of sets, ptrs and datasets
  op_set nodes(input->nnode, NULL);
  op_set edges(input->nedge, NULL);//&input->partition_length); //<THE CULPRIT!
  op_ptr edge_map(edges, nodes, 2, input->edge);
  op_dat<int> p_nodei(edges, 1, nodei);
  op_dat<int> p_nodej(edges, 1, nodej);
  op_dat<double> p_per_node(nodes, 3, per_node);
  op_dat<double> p_per_edge(edges, 3, per_edge);

  double *p_per_node_start[1];
  p_per_node_start[0] = (double *)p_per_node.dat_d;
  printf("(1) per_node_start[0] = %p\n", p_per_node_start[0]);
  op_dat_gbl<double *>p_p_per_node_start(1, p_per_node_start);

  //op_diagnostic_output();

  op_par_loop_5(res, edges,
    &p_nodei, 0, NULL, OP_READ,
    &p_nodej, 0, NULL, OP_READ,
    &p_p_per_node_start, 0, NULL, OP_READ,
    &p_per_node, 0, &edge_map, OP_INC,
    &p_per_edge, 0, NULL,      OP_RW);

  op_fetch_data((op_dat<double> *)&p_per_node);
  op_fetch_data((op_dat<double> *)&p_per_edge);

  // Check things have changed
  int num_per_node_zero = 0;
  for (int i=0; i<input->nnode*3; i++) {
    if (per_node[i] == 0.0) num_per_node_zero++;
  }
  int num_per_edge_zero = 0;
  for (int i=0; i<input->nedge*3; i++) {
    if (per_edge[i] == 0.0) num_per_edge_zero++;
  }
  printf("# num_per_node_zero %d/%d\n", num_per_node_zero, input->nnode*3);
  printf("# num_per_edge_zero %d/%d\n", num_per_edge_zero, input->nedge*3);

  //cleanup
  delete[] nodei;
  delete[] nodej;
  delete[] per_node;
  delete[] per_edge;
}
