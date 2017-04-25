function fxn_init() {
  // fxn_tribes(16)
  fxn_fitness_distrib_type_init()
  //fxn_selection_init();
  // assign values that were parsed from input file to JS vars
  // we need to store these, in case we need access them later
  tracking_threshold = dmi.tracking_threshold.value
  fraction_neutral = dmi.fraction_neutral.value
  // properly grey out items
  fxn_combine_mutns_able()
  fxn_dynamic_linkage_able()
  fxn_bottleneck_able()
  fxn_restart_case_able()
  // fxn_is_parallel_init()
  // fxn_migration()
  fxn_clone()
  fxn_fraction_neutral()
  //fxn_polygenic_beneficials();
  fxn_polygenics_able()
  fxn_track_neutrals()
  fxn_track_all_mutn()
  fxn_init_tracking_threshold()
  // show_hide_mutation_upload_form()
  fxn_auto_malloc()
  fxn_initial_alleles_init()
  //document.getElementById("tribediv").style.display = "none"
  //dmi.case_id.focus()
  compute_u()
}

function fxn_set_caseid() {
        dmi.case_id.value = dmi.cid.value
        //parent.frames.contents.caseidform.case_id.value = dmi.case_id.value
}

function fxn_initial_alleles() {
  if (dmi.initial_alleles.checked) {
    dmi.num_contrasting_alleles.readOnly = false
    dmi.max_total_fitness_increase.readOnly = false
    dmi.initial_alleles_pop_frac.readOnly = false
    dmi.initial_alleles_amp_factor.readOnly = false
    dmi.num_contrasting_alleles.value = "1000"
    dmi.max_total_fitness_increase.value = "1.0"
    dmi.initial_alleles_pop_frac.value = "1.0"
    dmi.num_contrasting_alleles.focus()
    $('#desc').tagsinput('add', 'Initial alleles');
  } else {
    dmi.num_contrasting_alleles.readOnly = true
    dmi.max_total_fitness_increase.readOnly = true
    dmi.initial_alleles_pop_frac.readOnly = true
    dmi.initial_alleles_amp_factor.readOnly = true
    dmi.num_contrasting_alleles.value = 0
    $('#desc').tagsinput('remove', 'Initial alleles');
  }
}

function fxn_initial_alleles_init() {
  if (dmi.num_contrasting_alleles.value == 0) {
    dmi.num_contrasting_alleles.readOnly = true
    dmi.max_total_fitness_increase.readOnly = true
    dmi.initial_alleles_pop_frac.readOnly = true
    dmi.initial_alleles_amp_factor.readOnly = true
  } else {
    dmi.num_contrasting_alleles.readOnly = false
    dmi.max_total_fitness_increase.readOnly = false
    dmi.initial_alleles_pop_frac.readOnly = false
    dmi.initial_alleles_amp_factor.readOnly = false
    dmi.initial_alleles.checked = true
  }
}

function validate(obj) {
  var val = parseFloat(obj.value)
  var max = parseFloat(obj.max)
  var min = parseFloat(obj.min)

  if (val < min || val > max) {
    obj.parentNode.parentNode.className = "form-group has-error"
  }else {
    obj.parentNode.parentNode.className = "form-group"
  }

}

function fxn_synergistic_epistasis() {
  fxn_synergistic_epistasis_able();
  if(dmi.synergistic_epistasis.checked) {
    if (dmi.se_nonlinked_scaling.value = "0.0"){
        dmi.se_nonlinked_scaling.value = "0.1"
    }
  }
}

function fxn_synergistic_epistasis_able() {
    if(dmi.synergistic_epistasis.checked) {
           dmi.se_nonlinked_scaling.readOnly = false
           dmi.se_linked_scaling.readOnly = false
    } else {
       dmi.se_nonlinked_scaling.readOnly = true
       dmi.se_linked_scaling.readOnly = true
    }
}
function fxn_synergistic_epistasis_disable() {
        dmi.se_nonlinked_scaling.readOnly = true
        dmi.se_linked_scaling.readOnly = true
}

function fxn_combine_mutns() {
  fxn_combine_mutns_able()
  if (dmi.combine_mutns.checked) {
    dmi.multiplicative_weighting.value = 0.5
    dmi.multiplicative_weighting.select()
  } else {
    dmi.multiplicative_weighting.value = 0.0
  }
}

function fxn_combine_mutns_able() {
  if (dmi.combine_mutns.checked) {
    document.getElementById("mwdiv").style.display = "block"
  } else {
    document.getElementById("mwdiv").style.display = "none"
  }
}

function fxn_dynamic_linkage() {
  fxn_dynamic_linkage_able()
  if (dmi.dynamic_linkage.checked) {
    //mendel_input.num_linkage_subunits.value = 989
    if (dmi.haploid_chromosome_number.value = "0")
           dmi.haploid_chromosome_number.value = "23"
  } else {
    dmi.num_linkage_subunits.value = 989
  }
}

function fxn_dynamic_linkage_able() {
    if (dmi.dynamic_linkage.checked) {
       dmi.haploid_chromosome_number.readOnly = false
           document.getElementById("num_linkage_subunits").innerText =
            "b. number of linkage subunits:"
    } else {
       dmi.haploid_chromosome_number.readOnly = true
       document.getElementById("num_linkage_subunits").innerText =
            "b. fixed block linkage number:"
    }
}

function fxn_haploid() {
   if (dmi.clonal_haploid.checked) {
      dmi.fraction_recessive.value = 0.0
      dmi.dominant_hetero_expression.value = 1.0
      status("Setting fraction_recessive to 0 and dominant_hetero_expression to 1")
   } else {
      dmi.dominant_hetero_expression.value = 0.5
      status("Setting dominant_hetero_expression back to 0.5")
   }
}

function fxn_is_parallel_init() {
   if (dmi.is_parallel.checked) {
      document.getElementById("psdiv").style.display = "block"
   } else {
      document.getElementById("psdiv").style.display = "none"
      document.getElementById("num_procs").value = 1
   }
}

function fxn_is_parallel() {
   tag = "Tribes"
   if (dmi.is_parallel.checked) {
      document.getElementById("psdiv").style.display = "block"
      np = document.getElementById("num_procs")
      if (np.value == 1) { np.value = 2 }
      $('#desc').tagsinput('add', tag);
   } else {
      document.getElementById("psdiv").style.display = "none"
      document.getElementById("num_procs").value = 1
      $('#desc').tagsinput('remove', tag);
   }
}

function status(msg) {
    try { // this will throw error if bootstrap-notify.min.js not included
        $.notify({ message: msg } , { placement: { from: "bottom", align: "right" }});
    } catch(err) { // the old way of doing things
        document.getElementById("warning").innerText = msg
    }
}

function warn(msg) {
    try {
        $.notify({ message: msg } , { type: 'warning', placement: { from: "bottom", align: "right" } });
    } catch(err) {
        document.getElementById("warning").innerText = msg
    }
}

function danger(msg) {
    try {
        $.notify({ message: msg } , { type: 'danger', placement: { from: "bottom", align: "right" } });
    } catch(err) {
        document.getElementById("danger").innerText = msg
    }
}

function fxn_restart_case() {
  fxn_restart_case_able()
  if (dmi.restart_case.checked) {
    if (dmi.restart_dump_number.value = "0") {
        dmi.restart_dump_number.value = "1"
    }
  }
}

function fxn_restart_case_able() {
   if (dmi.restart_case.checked) {
      document.getElementById("rddiv").style.display = "block"
   } else {
      document.getElementById("rddiv").style.display = "none"
   }
}

function fxn_bottleneck() {
  fxn_bottleneck_able()
  if (dmi.bottleneck_yes.checked) {
    if (dmi.bottleneck_generation.value = "0") {
        dmi.bottleneck_generation.value = "1000"
    }
    if (dmi.bottleneck_pop_size.value = "0") {
        dmi.bottleneck_pop_size.value = "100"
    }
    if (dmi.num_bottleneck_generations.value = "0") {
        dmi.num_bottleneck_generations.value = "500"
    }
  }
}
function fxn_bottleneck_able() {
   if (dmi.bottleneck_yes.checked) {
      document.getElementById("bydiv").style.display = "block"
      document.getElementById("nbg").style.display = "block"
   } else {
      document.getElementById("bydiv").style.display = "none"
   }
}

function check_bottleneck() {
   bgen = dmi.bottleneck_generation.value

   if(bgen < 0) {
     status("Cyclic bottlenecking turned on")
     if(dmi.num_bottleneck_generations.value > -bgen ) {
         dmi.num_bottleneck_generations.value = -bgen - 1;
     }
   } else {
     status("Cyclic bottlenecking turned off")
   }

   var pgr = dmi.pop_growth_rate

   e = document.getElementById("pop_growth_model")
   var selectedIndex = e.options[e.selectedIndex].value;

   if (selectedIndex == 4 && bgen < dmi.num_generations.value) {
      if (pgr.value - parseInt(pgr.value) == 0) {
        pgr.value = parseInt(pgr.value) + 0.2
        status("Updated pop growth rate to valid number")
      }
   }
}

function fxn_auto_malloc() {
  // estimate number of mutations that will be required
  // max number of mutations per individual
  compute_memory()

  var u = dmi.mutn_rate.value
  var uneu = u*dmi.fraction_neutral.value
  var uben = (u-uneu)*dmi.frac_fav_mutn.value
  var udel = (u-uneu)*(1-dmi.frac_fav_mutn.value)

  var ng = dmi.num_generations.value
  var min = 100
  var est_max_del = Math.max(ng*udel*2, min)
  var est_max_neu = Math.max(ng*uneu*2, min)
  var est_max_fav = Math.max(ng*uben*2, min)

  if (dmi.auto_malloc.checked) {
    dmi.max_del_mutn_per_indiv.value = est_max_del
    dmi.max_neu_mutn_per_indiv.value = est_max_neu
    dmi.max_fav_mutn_per_indiv.value = est_max_fav

    dmi.max_del_mutn_per_indiv.readOnly = true
    dmi.max_neu_mutn_per_indiv.readOnly = true
    dmi.max_fav_mutn_per_indiv.readOnly = true
  } else {
    dmi.max_del_mutn_per_indiv.readOnly = false
    dmi.max_neu_mutn_per_indiv.readOnly = false
    dmi.max_fav_mutn_per_indiv.readOnly = false
    dmi.max_del_mutn_per_indiv.select()
  }
}

function compute_memory() {
  // given in bytes
  var sizeof_float = 4;
  var sizeof_double = 8
  var sizeof_int = 4
  var sizeof_logical = 4

  // input variables affecting memory
  var opf = 2*dmi.reproductive_rate.value
  var frd = dmi.fraction_random_death.value
  var ng = dmi.num_generations.value
  var mutn_rate = dmi.mutn_rate.value
  var num_linkage_subunits = dmi.num_linkage_subunits.value
  var fraction_neutral = dmi.fraction_neutral.value
  var frac_fav_mutn = dmi.frac_fav_mutn.value
  if (dmi.pop_growth_model.value > 0) {
    pop_size = ng
  } else {
    pop_size = dmi.pop_size.value
  }

  var max_del = dmi.max_del_mutn_per_indiv.value
  var max_neu = dmi.max_neu_mutn_per_indiv.value
  var max_fav = dmi.max_fav_mutn_per_indiv.value

  // estimate memory requirements

  var u = dmi.mutn_rate.value
  var uneu = u*dmi.fraction_neutral.value
  var uben = (u-uneu)*dmi.frac_fav_mutn.value
  var udel = (u-uneu)*(1-dmi.frac_fav_mutn.value)

  var est_max_del = ng*udel
  var est_max_neu = ng*uneu
  var est_max_fav = ng*uben

  // highlight box as red if value too low

  if(!dmi.auto_malloc.checked) {
    if (max_del < est_max_del) {
      document.getElementById("max_del_mutn_per_indiv").className = "form-group has-error"
    } else {
      document.getElementById("max_del_mutn_per_indiv").className = "form-group"
    }
    if (max_neu < est_max_neu) {
      document.getElementById("max_neu_mutn_per_indiv").className = "form-group has-error"
    } else {
      document.getElementById("max_neu_mutn_per_indiv").className = "form-group"
    }
    if (max_fav < est_max_fav) {
      document.getElementById("max_fav_mutn_per_indiv").className = "form-group has-error"
    } else {
      document.getElementById("max_fav_mutn_per_indiv").className = "form-group"
    }
  }

  var max_size = 0.55*opf*(1.-frd)*pop_size
  var nmax = 12*(1. - frd) + 0.999

  //var base_mem = 174600; // bytes
  var base_mem = 174600; // bytes
  var mem_reqd = base_mem

  mem_reqd += num_linkage_subunits*6*max_size*sizeof_int;    // lb_mutn_count
  mem_reqd += num_linkage_subunits*2*max_size*sizeof_double; // lb_fitness
  mem_reqd += num_linkage_subunits*6*nmax*sizeof_int;    // offsprng_lb_mutn_count
  mem_reqd += num_linkage_subunits*2*nmax*sizeof_double; // offsprng_lb_fitness
  mem_reqd += max_size*sizeof_double;  // fitness
  mem_reqd += max_size*sizeof_double;  // fitness_score
  mem_reqd += max_size*sizeof_double;  // sorted_score
  mem_reqd += max_size*sizeof_logical; // available
  mem_reqd += pop_size*sizeof_logical; // replaced_by_offspring

  // for polymorphism analysis in diagnostics.f90
  mem_reqd += 1000000*sizeof_int

  // compute dynamic memory requirements
  dyn_mem_reqd =  max_del*max_size*sizeof_int; //dmutn
  dyn_mem_reqd += max_neu*max_size*sizeof_int; //nmutn
  dyn_mem_reqd += max_fav*max_size*sizeof_int; //fmutn
  dyn_mem_reqd += max_del*nmax*sizeof_int; //dmutn_offsprng
  dyn_mem_reqd += max_neu*nmax*sizeof_int; //nmutn_offsprng
  dyn_mem_reqd += max_fav*nmax*sizeof_int; //fmutn_offsprng

  // convert to MB
  mem_est = (dyn_mem_reqd + mem_reqd)/1024/1024
  var mem_est = Math.round(mem_est*100)/100

  msg = "Memory required for run: " + mem_est + " MB"
  document.getElementById("memory").innerText = msg

}

function compute_u() {
   u = dmi.mutn_rate.value
   uneu = u*dmi.fraction_neutral.value
   document.getElementById("uben").innerHTML = String(Math.round((u-uneu)*dmi.frac_fav_mutn.value*10)/10)
   document.getElementById("udel").innerHTML = String(Math.round((u-uneu)*(1-dmi.frac_fav_mutn.value)*10)/10)
   document.getElementById("uneu").innerHTML = String(Math.round(uneu*10)/10)
}

function fxn_fraction_neutral() {
   if(dmi.fraction_neutral.value > 0) {
      dmi.track_neutrals.checked = true
      //status("tracking all mutations")
   } else {
      dmi.track_neutrals.checked = false
      //status("not tracking all mutations")
   }
   compute_u()
}

function fxn_polygenics_able() {
  if(dmi.polygenic_beneficials.checked) {
    dmi.polygenic_init.readOnly = false
    dmi.polygenic_target.readOnly = false
    dmi.polygenic_effect.readOnly = false
  } else {
    dmi.polygenic_init.readOnly = true
    dmi.polygenic_target.readOnly = true
    dmi.polygenic_effect.readOnly = true
  }
}

function fxn_polygenic_beneficials(init) {
   fraction_neutral = dmi.fraction_neutral.value
   fraction_fav_mutn = dmi.frac_fav_mutn.value
   plot_allele_gens = dmi.plot_allele_gens.value
   tag = "Waiting time"
   if(dmi.polygenic_beneficials.checked) {
      dmi.polygenic_init.readOnly = false
      dmi.polygenic_target.readOnly = false
      dmi.polygenic_effect.readOnly = false
      dmi.track_neutrals.checked = true
      fxn_track_neutrals()
      dmi.fraction_neutral.value = 1.0
      dmi.frac_fav_mutn.value = 0.0
      dmi.dynamic_linkage.checked = false
      dmi.num_linkage_subunits.value = dmi.polygenic_target.value.length
      // document.getElementById("recombination_model").selectedIndex = 2
      document.getElementById("fitness_distrib_type").selectedIndex = 1
      fxn_fitness_distrib_type_init()
      dmi.uniform_fitness_effect_fav.readOnly = true
      dmi.haploid_chromosome_number.readOnly = true
      if(init==1) {
         dmi.plot_allele_gens.value = plot_allele_gens
      } else {
         dmi.plot_allele_gens.value = 1
      }
      compute_u()
      status("Turning on track_neutrals, setting fraction_neutral = 1.0, turning off dynamic linkage, setting num_linkage_subunits to length of target string, suppressing recombination, setting all mutations to equal effect")
      $('#desc').tagsinput('add', tag);
   } else {
      dmi.polygenic_init.readOnly = true
      dmi.polygenic_target.readOnly = true
      dmi.polygenic_effect.readOnly = true
      dmi.uniform_fitness_effect_del.readOnly = false
      dmi.haploid_chromosome_number.readOnly = true
      fxn_fitness_distrib_type_init()
      dmi.frac_fav_mutn.readOnly = false
      dmi.fraction_neutral.value = fraction_neutral
      dmi.frac_fav_mutn.value = fraction_fav_mutn
      dmi.plot_allele_gens.value = plot_allele_gens
      $('#desc').tagsinput('remove', tag);
   }
   fxn_auto_malloc()
}

function fxn_polygenic_target() {
   dmi.num_linkage_subunits.value = dmi.polygenic_target.value.length
   if(dmi.polygenic_init.value.length != dmi.polygenic_target.value.length) {
      danger("ERROR: polygenic init string must be same length as target")
      dmi.polygenic_init.select()
   }
}

function fxn_track_neutrals() {
   if(dmi.track_neutrals.checked) {
      dmi.fraction_neutral.readOnly = false
      if ( fraction_neutral > 0 ) {
         dmi.fraction_neutral.value = fraction_neutral
      } else {
         dmi.fraction_neutral.value = 0.9
      }
      // Modify mutation rate -- divide by fraction_neutrals
      //dmi.mutn_rate.value = Math.round(dmi.mutn_rate.value/(1-dmi.fraction_neutral.value))
      dmi.fraction_neutral.select()
      dmi.track_all_mutn.checked = true
      fxn_track_all_mutn()
      //document.getElementById("mutn_rate").innerText = "Total mutation rate per individual per generation:"
      warn("including neutrals in analysis will require more memory and will slow run, and all mutations will be tracked")
      $('#desc').tagsinput('add', 'Neutrals');
   } else {
      // Modify mutation rate -- multiply by fraction_neutrals
      //dmi.mutn_rate.value = Math.round(dmi.mutn_rate.value*(1-dmi.fraction_neutral.value))
      //document.getElementById("mutn_rate").innerText = "Total non-neutral mutation rate per individual per generation:"
      dmi.fraction_neutral.value = "0.0"
      dmi.fraction_neutral.readOnly = true
      $('#desc').tagsinput('remove', 'Neutrals');
   }
   compute_u()
}

function fxn_track_all_mutn() {
   if(dmi.track_all_mutn.checked) {
      dmi.tracking_threshold.value = 0
      dmi.tracking_threshold.readOnly = true
   } else {
      dmi.tracking_threshold.readOnly = false
      dmi.tracking_threshold.value = tracking_threshold
      dmi.tracking_threshold.select()
      dmi.track_neutrals.checked = false
   }
}

function fxn_init_tracking_threshold() {
   if(dmi.tracking_threshold.value == 0) {
      dmi.track_all_mutn.checked = true
   } else {
      dmi.track_all_mutn.checked = false
   }
   // set a default tracking_threshold value, in case user
   // is not restarting from a file that has a specified tracking
   // threshold value.
   if(tracking_threshold == 0) { tracking_threshold = 1.e-5; }
}

function fxn_fitness_distrib_type_init() {
   fdt = dmi.fitness_distrib_type.value
   // equal effect distribution
   if (fdt == 0) {
      document.getElementById("ufe_div").style.display = "block"
      document.getElementById("weibull_div").style.display = "none"
      document.getElementById("crdiv").style.display = "block"
      dmi.combine_mutns.readOnly = false
      dmi.synergistic_epistasis.readOnly = false
   // Weibull distribution
   } else if (fdt == 1) {
      document.getElementById("ufe_div").style.display = "none"
      document.getElementById("weibull_div").style.display = "block"
      document.getElementById("crdiv").style.display = "block"
      dmi.combine_mutns.readOnly = false
      dmi.synergistic_epistasis.readOnly = false
   // All Neutral
   } else if (fdt == 2) {
      document.getElementById("ufe_div").style.display = "none"
      document.getElementById("weibull_div").style.display = "none"
      document.getElementById("crdiv").style.display = "none"
      dmi.combine_mutns.readOnly = true
      dmi.synergistic_epistasis.readOnly = true
      fxn_disable_synergistic_epistasis()
   // Bi-modal
   } else if (fdt == 3) {
      document.getElementById("ufe_div").style.display = "block"
      document.getElementById("weibull_div").style.display = "block"
      document.getElementById("crdiv").style.display = "block"
      dmi.combine_mutns.readOnly = false
      dmi.synergistic_epistasis.readOnly = false
   } else {
      document.getElementById("ufe_div").style.display = "none"
      document.getElementById("weibull_div").style.display = "block"
      document.getElementById("crdiv").style.display = "block"
      dmi.combine_mutns.readOnly = false
      dmi.synergistic_epistasis.readOnly = false
   }
}

function fxn_fitness_distrib_type_change() {
   fxn_fitness_distrib_type_init()
   fdt = dmi.fitness_distrib_type.value
   if (fdt == 0) {
      dmi.dominant_hetero_expression.value = 1.0
      $('#desc').tagsinput('add', 'Equal effect mutations');
   } else {
      dmi.dominant_hetero_expression.value = 0.5
   }
}

function show_hide_advanced() {
    if (dmi.advsel.checked) {
           document.getElementById("tab-pane-1").style.display = "block"
        } else {
           document.getElementById("tab-pane-1").style.display = "none"
        }
}

function show_hide_parser() {
    if (document.parsed_data.show_data.checked) {
           document.getElementById("parser").style.display = "block"
        } else {
           document.getElementById("parser").style.display = "none"
        }
}

function show_hide_mutation_upload_form(i) {
    // if user checks upload mutations on the mutation pane
    // then automatically also check the upload mutations box
    // under population substructure, and vice-versa
    if(i==2) {
        if (dmi.altruistic.checked) {
            dmi.upload_mutations.checked = true
        } else {
            dmi.upload_mutations.checked = false
        }
    }

    // if user checks upload mutations on the mutation pane
    // then automatically also check the upload mutations box
    // under population substructure, and vice-versa
    if (dmi.upload_mutations.checked) {
        $("#upload_mutations_div").slideDown()
        $('#desc').tagsinput('add', 'Upload mutations');
        //dmi.mutn_file_id.readOnly = false
    } else if (dmi.altruistic.checked) {
        $("#upload_mutations_div").slideDown()
        $('#desc').tagsinput('add', 'Upload altruistic alleles');
    } else {
        $("#upload_mutations_div").slideUp()
        $('#desc').tagsinput('remove', 'Upload mutations');
        //dmi.mutn_file_id.readOnly = true
    }
}

function correct_lb() {
    var num_linkage_subunits = dmi.num_linkage_subunits.value
    var haploid_chromosome_number = dmi.haploid_chromosome_number.value
    var linkage_subunits_per_chromosome = parseInt(num_linkage_subunits / haploid_chromosome_number)
    // recompute num_linkage_subunits to be a integer multiplier of the number of chromosomes
    num_linkage_subunits = linkage_subunits_per_chromosome * haploid_chromosome_number
    dmi.num_linkage_subunits.value = num_linkage_subunits
}

function range(start, stop, step) {
  var a=[start], b=start;
  while(b < stop)  {b+=step; a.push(b)}
  return a;
}

function generate_mutations() {
    var pop_size = dmi.pop_size.value
    var num_linkage_subunits = dmi.num_linkage_subunits.value
    var haploid_chromosome_number = dmi.haploid_chromosome_number.value
    var linkage_subunits_per_chromosome = parseInt(num_linkage_subunits / haploid_chromosome_number)
    // recompute num_linkage_subunits to be a integer multiplier of the number of chromosomes
    num_linkage_subunits = linkage_subunits_per_chromosome * haploid_chromosome_number

    // handle combination of chromosome inputs, such as 1-3, 5-7, 21-23
    var chromosomes = $("#chromosome_range").val()
    console.log(chromosomes)

    nmutn = pop_size * num_linkage_subunits * 2
    if (nmutn > 10000) {
        $("#payload").text("Asking to generate " + nmutn + " mutations... too many.  Please reduce either pop_size or num_linkage_subunits such that pop_size x num_linkage_subunits x 2 <= 10000");
        return -1;
    }

    if (!chromosomes) {
        $("#payload").text("ERROR: must specify a chromosome number, e.g. 1")
        return -1;
    } else {
        $("#payload").text("")
    }

    var dominance = 1
    var frac_fav_mutn = dmi.frac_fav_mutn.value
    var genome_size = dmi.genome_size.value
    var max_fav_fitness_gain = dmi.max_fav_fitness_gain.value
    var high_impact_mutn_threshold = dmi.high_impact_mutn_threshold.value
    var high_impact_mutn_fraction = dmi.high_impact_mutn_fraction.value

    var alpha_del = Math.log(genome_size)
    var alpha_fav = Math.log(genome_size*max_fav_fitness_gain)
    var gamma_del = Math.log(-Math.log(high_impact_mutn_threshold)/alpha_del) /
                    Math.log(high_impact_mutn_fraction)
    var gamma_fav = Math.log(-Math.log(high_impact_mutn_threshold)/alpha_fav) /
                    Math.log(high_impact_mutn_fraction)
    // console.log(alpha_del)


    var fitness_distrib_type = dmi.fitness_distrib_type.value
    var chromo_array = chromosomes.split(',')
    console.log(chromo_array)

    new_chromo_array = Array()

    for (c in chromo_array) {
        x = chromo_array[c].split("-")
        console.log('x:' + x)
        start = parseInt(x[0])
        end = parseInt(x[1])
        if (end) {
            // new_chromo_array = range(start, end, 1)
            new_chromo_array.push.apply(new_chromo_array, range(start, end, 1))
        } else {
            new_chromo_array.push(start)
        }
    }
    console.log(new_chromo_array)

    var allele_id = 0
    for (i = 1; i <= pop_size; i++) {
        for (c in new_chromo_array) {
            chromosome = new_chromo_array[c]
            lb1 = (chromosome-1)*linkage_subunits_per_chromosome + 1
            lb2 = chromosome*linkage_subunits_per_chromosome
            for (lb = lb1; lb <= lb2; lb++) {
                x = Math.random()
                y = Math.random()
                if (y > frac_fav_mutn) {
                    if (fitness_distrib_type == 0) {
                        fitness = -dmi.uniform_fitness_effect_del.value*1.0
                    } else {
                        fitness = -Math.exp(-alpha_del*Math.pow(x, gamma_del))
                    }
                } else {
                    if (fitness_distrib_type == 0) {
                        fitness = dmi.uniform_fitness_effect_fav.value*1.0
                    } else {
                        fitness = max_fav_fitness_gain*Math.exp(-alpha_fav*Math.pow(x, gamma_fav))
                    }
                }
                for (haplotype = 1; haplotype <= 2; haplotype++) {
                    allele_id += 1
                    $("#payload").append(String(i) + ' ' + String(lb) + ' ' + String(haplotype) + ' ' + String(fitness.toExponential(8)) + ' ' + String(dominance)  + '&#xA;');
                }
            }
        }
    }

    if (allele_id > 3772) {
        $("#payload").text("ERROR: The number of alleles is too large to POST through the website.  Try reducing some combination of pop_size and num_linkage_subunits.  The current limit is set to 3772 alleles.");
        return -1;
    } else {
        $("#gen_stats").text("generated " + allele_id + " mutations")
    }

}

function fxn_migration() {
   x = parseInt(1*dmi.num_indiv_exchanged.value)
   max = parseInt(1*dmi.pop_size.value)
   if(x == 0) {
      dmi.migration_generations.readOnly = true
   } else {
      dmi.migration_generations.readOnly = false
      //dmi.migration_generations.value = 1
      if(x > max || x < 0) danger("ERROR: Value must be between 0 and " + max)
   }
   dmi.num_indiv_exchanged.value = x
}

function fxn_tribes(max_tribes) {
   myobject = dmi.num_procs
   num_procs = myobject.value

   // set max number of tribes for server from setting in config.inc
   if(num_procs > max_tribes) {
      myobject.value = max_tribes
      num_procs = max_tribes
   }
   // set min number of tribes
   //if(num_procs < 2) {
   //   myobject.value = 2;
   //   num_procs = 2
   //}

   //if (dmi.homogenous_tribes.checked) {
   //   document.getElementById("tribediv").style.display = "none"
   //} else {
   //   document.getElementById("tribediv").style.display = "block"
   //}

   if (dmi.tribal_competition.checked) {
      dmi.tc_scaling_factor.readOnly = false
      dmi.group_heritability.readOnly = false
      dmi.tc_scaling_factor.select()
      warn("Group competition is still under development. Proceed with caution. Must set extinction threshold > 0 for tribal fission.")
      if(dmi.extinction_threshold.value == 0.0)  {
         dmi.extinction_threshold.value = 0.1
      }
   } else {
      dmi.tc_scaling_factor.readOnly = true
      dmi.group_heritability.readOnly = true
   }

   dmi.num_procs.title = "2 - " + max_tribes
   // Add options to tribe_id select statement
   dmi.tribe_id.options.length=0
   for (i = 0; i < num_procs; i++) {
      a = (i+1)/1000 + ''; // compute number of tribe as a string
      b = a.substring(1);  // remove the leading 0 from 0.001, 0.002, etc.
      if((i+1)%10==0) b += '0'; // every 10 tribes: 0.01->0.010, 0.02->0.020, etc.
      dmi.tribe_id.options[i]=new Option(b, b, true, false)
   }
}

function fxn_clone() {
   if (dmi.recombination_model.selectedIndex == 2) {
      dmi.fraction_self_fertilization.readOnly = true
      dmi.num_contrasting_alleles.readOnly = true
      dmi.max_total_fitness_increase.readOnly = true
      dmi.dynamic_linkage.readOnly = true
      dmi.haploid_chromosome_number.readOnly = true
      dmi.num_linkage_subunits.value = 1
   } else {
      dmi.fraction_self_fertilization.readOnly = false
      dmi.num_contrasting_alleles.readOnly = false
      dmi.max_total_fitness_increase.readOnly = false
      dmi.dynamic_linkage.readOnly = false
      dmi.haploid_chromosome_number.readOnly = false
      //dmi.num_linkage_subunits.value = 1000
   }
}

function fxn_selection_init() {
  i = dmi.selection_scheme.value
  if (i == 1 || i == 2 || i == 4) {
     dmi.non_scaling_noise.value = 0.05
     //status("Setting non_scaling_noise to 0.05")
  } else {
     dmi.non_scaling_noise.value = 0.0
     //status("Setting non_scaling_noise to 0.0")
  }
}

function fxn_selection(i) {
  fxn_selection_init()
  if (i == 4) {
     document.getElementById("ptv").style.display = "block"
     dmi.partial_truncation_value.select()
  } else {
     document.getElementById("ptv").style.display = "none"
  }
}

function check_back_mutn() {
   tag = "Back mutations"
   if(dmi.allow_back_mutn.checked) {
      tt = dmi.tracking_threshold.value
      dmi.tracking_threshold.value = "0.0"
      status("Changed tracking threshold to 0.0 so that all mutations will be tracked")
      $('#desc').tagsinput('add', tag);
   } else {
      if(tt<=0) tt = 1.e-5;
      dmi.tracking_threshold.value = tt
      status("Changed tracking threshold back to " + tt )
      $('#desc').tagsinput('remove', tag);
   }
}

function fxn_pop_growth_model(i) {
  if (i == 0) {
     dmi.pop_growth_rate.readOnly = true
     dmi.carrying_capacity.readOnly = true
  } else if (i == 1) { // Exponential growth
     dmi.pop_growth_rate.readOnly = false
     dmi.carrying_capacity.readOnly = true
     dmi.pop_size.value = "2";
     dmi.num_generations.value = "2000";
     dmi.pop_growth_rate.value = "1.01";
     dmi.pop_growth_rate.title = "1.00 - 1.26";
     warn("WARNING: dynamic populations are experimental and largely untested")
     $('#desc').tagsinput('add', 'Exponential growth');
     $('#desc').tagsinput('remove', 'Carrying capacity');
     $('#desc').tagsinput('remove', 'Founder');
  } else if (i == 2) { // Carrying capacity
     dmi.pop_growth_rate.readOnly = false
     dmi.carrying_capacity.readOnly = false
     dmi.pop_size.value = "2";
     dmi.num_generations.value = "1000";
     dmi.pop_growth_rate.value = "0.1";
     dmi.pop_growth_rate.title = "0.0 - 1.0";
     warn("WARNING: dynamic populations are experimental and largely untested")
     $('#desc').tagsinput('add', 'Carrying capacity');
     $('#desc').tagsinput('remove', 'Exponential growth');
     $('#desc').tagsinput('remove', 'Founder');
  } else if (i == 3) { // Prescribed pop size
     dmi.pop_growth_rate.readOnly = true
     dmi.carrying_capacity.readOnly = true
     dmi.carrying_capacity.value = 10000
  } else if (i == 4) { // Founder effects
     dmi.pop_growth_rate.readOnly = false
     dmi.pop_growth_rate.value = "8"
     dmi.carrying_capacity.readOnly = false
     dmi.pop_size.value = "2";
     dmi.bottleneck_yes.checked = true
     dmi.bottleneck_generation.value = parseInt(dmi.num_generations.value) + 1
     dmi.bottleneck_pop_size.value = 10
     document.getElementById("bydiv").style.display = "block"
     document.getElementById("nbg").style.display = "none"
     $('#desc').tagsinput('add', 'Founder');
     $('#desc').tagsinput('remove', 'Carrying capacity');
     $('#desc').tagsinput('remove', 'Exponential growth');
  } else {
     dmi.pop_growth_rate.readOnly = false
  }
}
