#pragma once
#include <RcppThread.h>
#include <RcppEigen.h>
#include <algorithm>
#include <vinecopulib/misc/tools_stl.hpp>
#include <vinecopulib/bicop/class.hpp>

namespace vinecopulib {
  class Bicop;

  namespace tools_stl {
    template<typename T>
    std::vector<T> rank(const std::vector<T>& x)
    {
      std::vector<size_t> r(x.size());
      auto order = tools_stl::get_order(x);
      for (auto i : order)
        r[order[i]] = static_cast<double>(i + 1);
      return r;
    }
  }
}

namespace vinereg {
using namespace vinecopulib;

struct CVineFitTemporaries
{
  std::vector<Eigen::VectorXd> hfunc1a;
  std::vector<Eigen::VectorXd> hfunc1b;
  std::vector<Eigen::VectorXd> hfunc1a_sub;
  std::vector<Eigen::VectorXd> hfunc1b_sub;
  std::vector<Bicop> pcs;
  std::vector<size_t> remaining_vars;
  std::vector<size_t> selected_vars;
  double crit;
};

class CVineRegSelector
{
public:
  CVineRegSelector(const Eigen::MatrixXd& data,
                   const std::vector<std::string>& var_types,
                   const FitControlsBicop& controls);

  void select_model();
  std::vector<size_t> get_selected_vars() const { return fit_.selected_vars; }
  std::vector<std::vector<Bicop>> get_pcs() const { return pcs_; }

private:
  void extend_fit(CVineFitTemporaries& fit, size_t var) const;
  void initialize_var(CVineFitTemporaries& fit, size_t var) const;
  std::vector<std::string> get_edge_types(const CVineFitTemporaries& fit,
                                          size_t t) const;
  Eigen::MatrixXd get_edge_data(const CVineFitTemporaries& fit, size_t t) const;
  void fit_pair_copula(CVineFitTemporaries& fit,
                       size_t t,
                       const Eigen::MatrixXd& u_e) const;
  void update_hfunc1a(CVineFitTemporaries& fit,
                     size_t t,
                     const Eigen::MatrixXd& u_e) const;
  void update_hfunc1b(CVineFitTemporaries& fit,
                     size_t t,
                     const Eigen::MatrixXd& u_e) const;
  void update_hfuncs(CVineFitTemporaries& fit,
                     size_t t,
                     const Eigen::MatrixXd& u_e) const;
  void update_selcrit(CVineFitTemporaries& fit) const;
  void update_vars(CVineFitTemporaries& fit, size_t var) const;
  void update_status(CVineFitTemporaries& fit, size_t var) const;

  size_t p_;
  Eigen::MatrixXd data_;
  std::vector<std::string> var_types_;
  FitControlsBicop controls_;
  CVineFitTemporaries fit_;
  std::vector<std::vector<Bicop>> pcs_;
};

inline CVineRegSelector::CVineRegSelector(
    const Eigen::MatrixXd& data,
    const std::vector<std::string>& var_types,
    const FitControlsBicop& controls)
  : p_(var_types.size() - 1)
  , data_(data)
  , var_types_(var_types)
  , controls_(controls)
{
  fit_.hfunc1a.resize(p_ + 1); // p+1, as otherwise error in update_hfuncs
  fit_.hfunc1b.resize(p_ + 1); // p+1, as otherwise error in update_hfuncs
  fit_.hfunc1a_sub.resize(p_); // needs discrete adaption
  fit_.hfunc1b_sub.resize(p_); // needs discrete adaption
  fit_.pcs.resize(p_);
  fit_.remaining_vars = tools_stl::seq_int(1, p_);
  fit_.selected_vars.reserve(p_);
  fit_.crit = 0.0;

  fit_.hfunc1a[0] = data_.col(0);

  // needs discrete adaption

}

inline void CVineRegSelector::select_model()
{
  std::mutex m; // required to synchronize write/reads to the selector
  auto num_threads = controls_.get_num_threads();
  RcppThread::ThreadPool pool(num_threads > 1 ? num_threads : 0);
  controls_.set_num_threads(0);

  while (fit_.selected_vars.size() < p_) {
    auto old_fit = fit_;  // fix current model (fit_ will be modified below)
    auto fit_replace_if_better = [&](size_t var) {
      CVineFitTemporaries new_fit = old_fit;
      this->extend_fit(new_fit, var);
      std::lock_guard<std::mutex> lk(m);  // synchronize
      if (new_fit.crit > fit_.crit)
        fit_ = std::move(new_fit);
    };
    pool.map(fit_replace_if_better, old_fit.remaining_vars);
    pool.wait();

    if (fit_.selected_vars == old_fit.selected_vars)
      break;  // could not improve the selection criterion

    // model improved; store pair copulas for new variable
    auto p_sel = fit_.selected_vars.size();
    pcs_.push_back(std::vector<Bicop>{ fit_.pcs[p_sel - 1] });
    for (size_t t = 0; t < p_sel - 1; t++)
      pcs_[t].push_back(fit_.pcs[t]);
  }
  pool.join();
  controls_.set_num_threads(num_threads);
}

inline void CVineRegSelector::extend_fit(CVineFitTemporaries& fit,
                                         size_t var) const
{
  this->initialize_var(fit, var);

  for (size_t t = 0; t < fit.selected_vars.size() + 1; t++) {
    auto u_e = this->get_edge_data(fit, t);
    this->fit_pair_copula(fit, t, u_e);
    this->update_hfuncs(fit, t, u_e);
  }
  this->update_status(fit, var);
}

inline void CVineRegSelector::initialize_var(CVineFitTemporaries& fit,
                                             size_t var) const
{
  fit.hfunc1b[0] = data_.col(var);
  // needs discrete adaption
}

// obtain variable types for the new edge in tree t
inline std::vector<std::string> CVineRegSelector::get_edge_types(
    const CVineFitTemporaries& fit, size_t t) const
{
  // the variable type can be inferred from the existence of _sub data
  std::vector<std::string> var_types(2);
  var_types[0] = "c";  // needs discrete adaption
  var_types[1] = "c"; // needs discrete adaption
  return var_types;
}

// obtain data for the new edge in tree t
inline Eigen::MatrixXd CVineRegSelector::get_edge_data(
    const CVineFitTemporaries& fit, size_t t) const
{
  Eigen::MatrixXd u_e(data_.rows(), 2);

  // not all trees have been fit
  if (t != fit.selected_vars.size()) {
    u_e.col(0) = fit.hfunc1a[t];
    u_e.col(1) = fit.hfunc1b[t];
  } else {  // last tree fit
    u_e.col(0) = fit.hfunc1b[t];
    u_e.col(1) = fit.hfunc1a[t];
  }

  // needs discrete adaption

  return u_e;
}

inline void CVineRegSelector::fit_pair_copula(CVineFitTemporaries& fit,
                                              size_t t,
                                              const Eigen::MatrixXd& u_e) const
{
  auto var_types = this->get_edge_types(fit, t);
  fit.pcs[t].set_var_types(var_types);
  fit.pcs[t].select(u_e, controls_);
}

inline void CVineRegSelector::update_hfunc1a(CVineFitTemporaries& fit,
                                            size_t t,
                                            const Eigen::MatrixXd& u_e) const
{

  // last tree fit
  if (t == fit.selected_vars.size()) {
    fit.hfunc1a[t] = fit.hfunc1b[t];
    fit.hfunc1a[t + 1] = fit.pcs[t].hfunc1(u_e);
  }

  // needs discrete adaption

}

inline void CVineRegSelector::update_hfunc1b(CVineFitTemporaries& fit,
                                            size_t t,
                                            const Eigen::MatrixXd& u_e) const
{

  // not all trees have been fit
  if (t != fit.selected_vars.size()) {
    fit.hfunc1b[t + 1] = fit.pcs[t].hfunc1(u_e);
  }

  // needs discrete adaption

}

void CVineRegSelector::update_hfuncs(CVineFitTemporaries& fit,
                                     size_t t,
                                     const Eigen::MatrixXd& u_e) const
{
  this->update_hfunc1a(fit, t, u_e);
  this->update_hfunc1b(fit, t, u_e);
}


// update value of the criterion used for variable and family selection
inline void CVineRegSelector::update_selcrit(CVineFitTemporaries& fit) const
{
  if (controls_.get_selection_criterion() == "loglik")
    fit.crit += fit.pcs[fit.selected_vars.size()].get_loglik();
  if (controls_.get_selection_criterion() == "aic")
    fit.crit -= fit.pcs[fit.selected_vars.size()].get_aic();
  if (controls_.get_selection_criterion() == "bic")
    fit.crit -= fit.pcs[fit.selected_vars.size()].get_bic();
}

// remove var from remaining variables; add to selected variables
inline void CVineRegSelector::update_vars(CVineFitTemporaries& fit, size_t var)
  const
{
  fit.remaining_vars.erase(
    std::remove(fit.remaining_vars.begin(), fit.remaining_vars.end(), var));
  fit.selected_vars.push_back(var);
}

inline void CVineRegSelector::update_status(CVineFitTemporaries& fit,size_t var)
  const
{
  this->update_selcrit(fit);
  this->update_vars(fit, var);
}

}
