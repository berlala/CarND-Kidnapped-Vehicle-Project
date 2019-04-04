/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 * Update April 4, 2019 by Bolin
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
using std::uniform_real_distribution;
using std::uniform_int_distribution;


void ParticleFilter::init(double x, double y, double theta, double std_pos[]) {
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
  particles.resize(num_particles);  // add after first review

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
			double x_predicted = predicted[ii].x;
			double y_predicted = predicted[ii].y;
			double dist_c = dist(x_obs, y_obs, x_predicted, y_predicted);
			//double dist_single_point[ii] = { dist_c };
			dist_single_point.push_back(dist_c);//modify after review
		}
		auto dist_smallest = std::min_element(std::begin(dist_single_point), std::end(dist_single_point));
		int index_smallest = std::distance(std::begin(dist_single_point), dist_smallest);
		//observations[i].id = index_smallest; // 该观测对应的index值
		observations[i].id = predicted[index_smallest].id;//modify after review
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

	for (int i = 0; i < num_particles; i++) {

		double x_p = particles[i].x;
		double y_p = particles[i].y;
		double theta_p = particles[i].theta;

		// create a vector to hold the map landmark locations predicted to be within sensor range of the particle
		vector<LandmarkObs> pred_map;

		// for each map landmark...
		for (unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++) {

			// get id and x,y coordinates
			float x_lm_map = map_landmarks.landmark_list[j].x_f;
			float y_lm_map = map_landmarks.landmark_list[j].y_f;
			int id_lm_map = map_landmarks.landmark_list[j].id_i;

			// only consider landmarks within sensor range of the particle (rather than using the "dist" method considering a circular 
			// region around the particle, this considers a rectangular region but is computationally faster)
			if (fabs(x_p - x_lm_map) <= sensor_range && fabs(y_p - y_lm_map) <= sensor_range) {

				// add prediction to vector
				pred_map.push_back(LandmarkObs{ id_lm_map, x_lm_map, y_lm_map }); // the effective landmarks in the sensor range
			}
		}

		// create and populate a copy of the list of observations transformed from vehicle coordinates to map coordinates
		vector<LandmarkObs> obs_map;
		for (unsigned int j = 0; j < observations.size(); j++) {
			double x_t = cos(theta_p)*observations[j].x - sin(theta_p)*observations[j].y + x_p;
			double y_t = sin(theta_p)*observations[j].x + cos(theta_p)*observations[j].y + y_p;
			obs_map.push_back(LandmarkObs{ observations[j].id, x_t, y_t });
		}

		// perform dataAssociation for the predictions and transformed observations on current particle
		dataAssociation(pred_map, obs_map);//关联《应该被观测到的landmark在地图上的位置》和《实际基于该粒子位置加观测得到的landmark在地图上的位置》

		// init weight
		particles[i].weight = 1.0;

		for (unsigned int j = 0; j < obs_map.size(); j++) {

			double x_o, y_o, x_pred_ass, y_pred_ass;
			x_o = obs_map[j].x;
			y_o = obs_map[j].y;
			int ass_id = obs_map[j].id;

			// get the x,y coordinates of the prediction associated with the current observation
			for (unsigned int k = 0; k < pred_map.size(); k++) {
				if (pred_map[k].id == ass_id) {
					x_pred_ass = pred_map[k].x;
					y_pred_ass = pred_map[k].y;
				}
			}

			//多维高斯分布计算
			double x_std = std_landmark[0];
			double y_std = std_landmark[1];
			double weight_obs = (1 / (2 * M_PI*x_std*y_std)) * exp(-(pow(x_pred_ass - x_o, 2) / (2 * pow(x_std, 2)) + (pow(y_pred_ass - y_o, 2) / (2 * pow(y_std, 2)))));

			// product of this obersvation weight with total observations weight
			particles[i].weight *= weight_obs;
		}
	} // end for each particle
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
	vector<Particle> new_particles;

	// get all of the current weights
	vector<double> weights;
	for (int i = 0; i < num_particles; i++) {
		weights.push_back(particles[i].weight);
	}

	//生成随机index
	default_random_engine rnd;
	uniform_int_distribution<int> rand_index(0, num_particles - 1);
	auto index = rand_index(rnd);

	// 最大权重指针iterator pointing
	double max_weight = *max_element(weights.begin(), weights.end());

	// uniform random distribution [0.0, max_weight)
	uniform_real_distribution<double> rand_weight(0.0, max_weight);

	double beta = 0.0;

	//更新粒子群
	for (int i = 0; i < num_particles; i++) {
		beta += rand_weight(rnd) * 2.0;
		while (beta > weights[index]) {
			beta -= weights[index];
			index = (index + 1) % num_particles;
		}
		new_particles.push_back(particles[index]);
	}

	particles = new_particles;
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