/*
 * OP2 instantiation of the hertz pairwise kernel.
 *
 * We take two datasets (from file) inputs (per-particle and contact data) and
 * expected output following a pairwise calculation.
 *
 * We use OP2 to express the problem as an iteration over contacts.
 *
 */

//#define AOS_EXTRA_DEBUG //< add (i,j) index information to struct

#ifdef GPU_TIMER
  #include "cuda_timer.h"
#else
  #include "simple_timer.h"
#endif
#include "check_result_vector.h"
#include "hertz_constants.h"
#include "user_defined_types.h"
#include "op_seq.h"
#include "res.h"
#include "update.h"
#include "framework.h"
#include <sstream>

using namespace std;

SimpleTimer dummy[8];

void run(struct params *input, int num_iter) {

  //--------------------
  // One-time only costs
  //--------------------

  //AoS generate
  one_time.push_back(SimpleTimer("aos_gen"));
  one_time.back().start();
  struct contact *aos = new contact[input->nedge];
  for (int e=0; e<input->nedge; e++) {
    int i = input->edge[(e*2)];
    int j = input->edge[(e*2)+1];

#ifdef AOS_EXTRA_DEBUG
    aos[e].i = i;
    aos[e].j = j;
#endif

    aos[e].xi[0] = input->x[(i*3)];
    aos[e].xi[1] = input->x[(i*3)+1];
    aos[e].xi[2] = input->x[(i*3)+2];
    aos[e].xj[0] = input->x[(j*3)];
    aos[e].xj[1] = input->x[(j*3)+1];
    aos[e].xj[2] = input->x[(j*3)+2];

    aos[e].vi[0] = input->v[(i*3)];
    aos[e].vi[1] = input->v[(i*3)+1];
    aos[e].vi[2] = input->v[(i*3)+2];
    aos[e].vj[0] = input->v[(j*3)];
    aos[e].vj[1] = input->v[(j*3)+1];
    aos[e].vj[2] = input->v[(j*3)+2];

    aos[e].omegai[0] = input->omega[(i*3)];
    aos[e].omegai[1] = input->omega[(i*3)+1];
    aos[e].omegai[2] = input->omega[(i*3)+2];
    aos[e].omegaj[0] = input->omega[(j*3)];
    aos[e].omegaj[1] = input->omega[(j*3)+1];
    aos[e].omegaj[2] = input->omega[(j*3)+2];

    aos[e].radiusi = input->radius[i];
    aos[e].radiusj = input->radius[j];

    aos[e].massi = input->mass[i];
    aos[e].massj = input->mass[j];

    aos[e].typei = input->type[i];
    aos[e].typej = input->type[j];
  }
  double *force_inc = new double[input->nnode*3];
  double *force_dec = new double[input->nnode*3];
  double *torque_dec = new double[input->nnode*3];
  one_time.back().stop_and_add_to_total();

  // OP2 init
  const char *fake_argv = {"foo"};
  one_time.push_back(SimpleTimer("op2_init"));
  one_time.back().start();
  op_init(1, (char **)fake_argv);
  one_time.back().stop_and_add_to_total();

  // OP2 declarations of sets, ptrs and datasets
  one_time.push_back(SimpleTimer("op2_decl"));
  one_time.back().start();
  op_set nodes(input->nnode, NULL);
  op_set edges(input->nedge, (input->npartition == 1 ? NULL : &input->partition_length));
  op_ptr edge_map(edges, nodes, 2, input->edge);
  op_dat<struct contact> p_aos(edges, 1, aos);
  op_dat<double> p_force_inc(nodes, 3, force_inc);
  op_dat<double> p_force_dec(nodes, 3, force_dec);
  op_dat<double> p_torque_dec(nodes, 3, torque_dec);
  op_dat<double> p_shear(edges, 3, input->shear);
  op_dat<double> p_force(edges, 3, input->force);
  op_dat<double> p_torque(edges, 3, input->torque);
  one_time.back().stop_and_add_to_total();

  // OP2 plan (see HACK below)
  one_time.push_back(SimpleTimer("op2_plan"));

  op_diagnostic_output();

  //--------------------
  // Per-iteration costs
  //--------------------
  double res_time[2];
  per_iter.push_back(SimpleTimer("compute_kernel"));
  per_iter.push_back(SimpleTimer("add_kernel"));
  per_iter.push_back(SimpleTimer("result_fetch"));
  for (int run=0; run<num_iter; run++) {

    //PREPROCESSING
    // prep for run (not necessary to time these)
    for (int i=0; i<input->nnode*3; i++) {
      force_inc[i] = 0.0;
      force_dec[i] = 0.0;
      torque_dec[i] = 0.0;
    }

    //-----------------------------------------------------------------------

    //KERNEL INVOCATION
    //res kernel computes delta values (force/torque) and shear
    //update kernel produces force/torque results
    per_iter[0].start();
    op_par_loop_6(res, edges,
        &p_aos,    0, NULL, OP_READ,
        &p_shear,  0, NULL, OP_RW,
        &p_force_inc, 0, &edge_map, OP_INC,
        &p_force_dec, 1, &edge_map, OP_INC,
        &p_torque_dec, 0, &edge_map, OP_INC,
        &p_torque_dec, 1, &edge_map, OP_INC);
    per_iter[0].stop_and_add_to_total();
    //HACK for op2_plan one-time cost
    //This cost is wrapped up into the first invocation of the res kernel 
    //(run 0). We record the runtime for the first two invocations of the
    //kernel and use the difference to approximate the plan cost. We also fixup
    //the per_iter cost for the res kernel.
    if (run == 0) {
      res_time[0] = per_iter[0].total_time();
    } else if (run == 1) {
      //total_time is a cumulative total so we subtract the run 0 time
      res_time[1] = per_iter[0].total_time() - res_time[0];
      double op2_plan_time = res_time[0] - res_time[1];
      one_time.back().set_total_time(op2_plan_time);
      per_iter[0].set_total_time(2*res_time[1]);
    }

    op_fetch_data((op_dat<double> *)&p_shear);
    op_fetch_data((op_dat<double> *)&p_force_inc);
    op_fetch_data((op_dat<double> *)&p_force_dec);
    op_fetch_data((op_dat<double> *)&p_torque_dec);

    //emit first 3 values of force_inc/force_dec
    for (int i=0; i<9; i++)
      printf("%sforce_inc[%d] = %.16f\n", (i % 3 == 0) ? "\n" : "", i, force_inc[i]);
    for (int i=0; i<9; i++)
      printf("%sforce_dec[%d] = %.16f\n", (i % 3 == 0) ? "\n" : "", i, force_dec[i]);
    for (int i=0; i<9; i++)
      printf("%sinput->force[%d] = %.16f\n", (i % 3 == 0) ? "\n" : "", i, input->force[i]);

/*
    per_iter[1].start();
    op_par_loop_5(update, nodes,
        &p_force_inc, 0, NULL, OP_READ,
        &p_force_dec, 0, NULL, OP_READ,
        &p_torque_dec, 0, NULL, OP_READ,
        &p_force, 0, NULL, OP_RW,
        &p_torque, 0, NULL, OP_RW);
    per_iter[1].stop_and_add_to_total();

    //-----------------------------------------------------------------------

    //POSTPROCESSING
    //memcpy data back to host
    per_iter[2].start();
    op_fetch_data((op_dat<double> *)&p_shear);
    op_fetch_data((op_dat<double> *)&p_force);
    op_fetch_data((op_dat<double> *)&p_torque);
    per_iter[2].stop_and_add_to_total();
*/

    //-----------------------------------------------------------------------

#if 1
    //CHECKING
    //only check results the first time around
    if (run == 0) {
      for (int n=0; n<input->nnode; n++) {
        std::stringstream out;
        out << "force[" << n << "]";
        check_result_vector(
            out.str().c_str(),
            &input->expected_force[(n*3)], &input->force[(n*3)], 0.00001, true);
        out.str("");

      //out << "torque[" << n << "]";
      //check_result_vector(
      //    out.str().c_str(),
      //    &input->expected_torque[(n*3)], &input->torque[(n*3)]);
      }
      for (int n=0; n<input->nedge; n++) {
        stringstream out;
        out << "shear[" << n << "]";
        check_result_vector(
            out.str().c_str(), 
            &input->expected_shear[(n*3)], &input->shear[(n*3)]);
      }
    }
#endif

  }

  //cleanup
  delete[] aos;
  delete[] force_inc;
  delete[] force_dec;
  delete[] torque_dec;
}
