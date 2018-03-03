/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

    default_random_engine gen;

    // Lines below create a normal (Gaussian) distribution for x,y & theta.
    normal_distribution<double> dist_x(0, std[0]);
    normal_distribution<double> dist_y(0, std[1]);
    normal_distribution<double> dist_theta(0, std[2]);

    //Set number of particles
    num_particles = 100;

    Particle part;
    for(int i=0 ; i < num_particles; i++){
        part.id = i;
        part.x = x + dist_x(gen);
        part.y = y + dist_y(gen);
        part.theta = theta + dist_theta(gen);
        part.weight = 1.0;
        particles.push_back(part);
    }

    is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	default_random_engine gen;

    normal_distribution<double> dist_x(0, std_pos[0]);
    normal_distribution<double> dist_y(0, std_pos[1]);
    normal_distribution<double> dist_theta(0, std_pos[2]);


    for (int i = 0 ; i < num_particles; i++){
        if(fabs(yaw_rate) < 0.00001){
            particles[i].x+= cos(particles[i].theta)*delta_t*velocity;
            particles[i].y+= sin(particles[i].theta)*delta_t*velocity;
        }
        else{
            particles[i].x+= velocity/yaw_rate*(sin(particles[i].theta+yaw_rate*delta_t)-sin(particles[i].theta));
            particles[i].y+= velocity/yaw_rate*(cos(particles[i].theta)-cos(particles[i].theta+yaw_rate*delta_t));
            particles[i].theta+= yaw_rate*delta_t;
        }
        particles[i].x+= dist_x(gen);
        particles[i].y+= dist_y(gen);
        particles[i].theta+= dist_theta(gen);
    }


}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

    for(int i= 0; i <observations.size(); i++){

        int mapId = -1;

        double distance = numeric_limits<double>::max();

        for(int j = 0; j<predicted.size(); j++){
            double closest = dist(observations[i].x,observations[i].y,predicted[j].x,predicted[j].y);
            if(closest< distance){
                distance = closest;
                mapId = predicted[j].id;
            }
        }

        observations[i].id = mapId;
    }


}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

    for(int i = 0 ; i < particles.size(); i++){
        // Vector of predicted measurements
        std::vector<LandmarkObs> predMnt;
        for(int j = 0 ; j< map_landmarks.landmark_list.size(); j++){
            // check if the predicted measurements to all map within the sensor range for each particle
            if(fabs(particles[i].x-map_landmarks.landmark_list[j].x_f)<=sensor_range &&
               fabs(particles[i].y-map_landmarks.landmark_list[j].y_f)<=sensor_range) {
                //placeholder structure
                LandmarkObs p;
                p.id = map_landmarks.landmark_list[j].id_i;
                p.x = map_landmarks.landmark_list[j].x_f;
                p.y = map_landmarks.landmark_list[j].y_f;
                predMnt.push_back(p);
            }
        }


        //transform observations
        std::vector<LandmarkObs> trnsfObvs;
        particles[i].weight=1.0;

        for(int k = 0; k < observations.size(); k++){
            LandmarkObs obs;
            obs.x = particles[i].x+(cos(particles[i].theta)*observations[k].x) -
                    (sin(particles[i].theta)*observations[k].y);
            obs.y = particles[i].y+(sin(particles[i].theta)*observations[k].x) +
                    (cos(particles[i].theta)*observations[k].y);
            obs.id = observations[k].id;
            trnsfObvs.push_back(obs);
        }

        dataAssociation(predMnt,trnsfObvs);
        double x_obs, y_obs, mu_x, mu_y;
        for (int k = 0; k <trnsfObvs.size(); k++){
            x_obs = trnsfObvs[k].x;
            y_obs = trnsfObvs[k].y;
            for(int l = 0 ; l <predMnt.size(); l++){
                if(trnsfObvs[k].id == predMnt[l].id){
                    mu_x = predMnt[l].x;
                    mu_y = predMnt[l].y;
                }
            }
            double gauss_norm = (1/(2*M_PI*std_landmark[0]*std_landmark[1]));
            double exponent = (pow((x_obs - mu_x),2)/(2*pow(std_landmark[0],2))) +
                              (pow((y_obs - mu_y),2)/(2*pow(std_landmark[1],2)));
            double particleWt = gauss_norm * exp(-exponent);
            particles[i].weight*= particleWt;
        }

    }

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

    default_random_engine gen;

    vector<Particle> resampled;

    vector<double>wtOnly;

    for(int i=0; i<num_particles; i++){
        wtOnly.push_back(particles[i].weight);
    }

    uniform_int_distribution<int> index(0,num_particles-1);
    unsigned int current = index(gen);

    double beta = 0.0;

    double maxWt = *max_element(wtOnly.begin(), wtOnly.end());
    uniform_real_distribution<double> randWt(0.0, maxWt);

    for (int i = 0 ; i < num_particles; i ++){
        beta += 2.0*randWt(gen);

        while (beta > wtOnly[current]) {
            beta -= wtOnly[current];
            current = (current + 1) % num_particles;
        }
        resampled.push_back(particles[current]);
    }
    particles = resampled;

}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
