static void update(
  double *force_inc, double *force_dec,
  double *torque_dec,
  double *force, 
  double *torque) {
  force[0] += force_inc[0];
  force[1] += force_inc[1];
  force[2] += force_inc[2];

  force[0] -= force_dec[0];
  force[1] -= force_dec[1];
  force[2] -= force_dec[2];

  torque[0] -= torque_dec[0];
  torque[1] -= torque_dec[1];
  torque[2] -= torque_dec[2];
}
