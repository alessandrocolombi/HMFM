//
// Created by pietr on 12/12/2021.
//

#include "GibbsSampler.h"


void GibbsSampler::GS_Step() {
    for(FullConditional* full_cond: this->FullConditionals){ // mettere prima update della partition (da aggiungere anche prima)
        full_cond->update(gs_data, random_engine);
    }
}


void GibbsSampler::store_params_values() {
    for(std::map<string,std::vector<float>>::iterator iter = output_data.begin(); iter != output_data.end(); ++iter)
        output_data[iter->first].append(parameters[iter->first]);
}


std::map<string, std::vector<float>> GibbsSampler::sample() {
    for(unsigned int it; it<burn_in + n_iter * thin; it++){
        this->GS_Step();
        if(it>burn_in && it%thin == 0){
            this->store_params_values();
        }
    }
    return output_data;
}