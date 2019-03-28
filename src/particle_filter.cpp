/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;
using std::default_random_engine;


void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 100;  // TODO: Set the number of particles
  //add by bolin 2019/3/19
  default_random_engine rnd;
  normal_distribution<double> x_init(x, std_pos[0]);
  normal_distribution<double> y_init(y, std_pos[1]);
  normal_distribution<double> theta_init(theta, std_pos[2]);

  for (int a = 0; a < num_particles; ++a) {
	  particles[a].x = x_init(rnd);
	  particles[a].y = y_init(rnd);
	  particles[a].theta = theta_init(rnd);
	  particles[a].id = a;
	  particles[a].weight = 1;
  }
  bool init_flag = 1;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
	default_random_engine rnd;
	for (int i = 0; i < num_particles; ++i) {
		double x_pred;
		double y_pred;
		double theta_pred;
		if (yaw_rate < 0.0001) {
			x_pred = particles[i].x + velocity * cos(particles[i].theta)*delta_t;
			y_pred = particles[i].y + velocity * sin(particles[i].theta)*delta_t;
			theta_pred = particles[i].theta;
		}
		else {
			x_pred= particles[i].x + (velocity / yaw_rate) * (sin(particles[i].theta + (yaw_rate * delta_t)) - sin(particles[i].theta));
			y_pred = particles[i].y + (velocity / yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta + (yaw_rate * delta_t)));
			theta_pred = particles[i].theta + (yaw_rate * delta_t); //ref to Lesson 6.8
		}
		normal_distribution<double> x_gpred(x_pred, std_pos[0]);
		normal_distribution<double> y_gpred(y_pred, std_pos[1]);
		normal_distribution<double> theta_gpred(theta_pred, std_pos[2]);

		particles[i].x = x_gpred(rnd);
		particles[i].y = y_gpred(rnd);
		particles[i].theta = theta_gpred(rnd);
	}

	bool predict_flag = 1;

}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
	for (int i = 0; i < observations.size(); ++i) { //对每一个观测值
		double x_obs = observations[i].x;
		double y_obs = observations[i].y;
		std::vector<double> dist_single_point;
		for (int ii = 0; ii < predicted.size(); ++ii) { // 对每一个predict的值（为什么是predict）
			double dist_c = dist(x_obs, y_obs, double predicted[ii].x, double predicted[ii].y)；
			double dist_single_point[ii] = dist_c;
		}
		auto dist_smallest = std::min_element(std::begin(dist_single_point), std::end(dist_single_point));
		int index_smallest = std::distance(std::begin(dist_single_point), dist_single_point);
		int observations[i].id = index_smallest; // 该观测对应的index值
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}