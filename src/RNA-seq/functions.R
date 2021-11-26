## ---- create_functions
sample_contrib_plot <- function(df,x,y,facet_factor,rows){
  ggplot(data=df, aes(x=.data[[x]], y=.data[[y]], fill = .data[[x]]))+
    geom_bar(stat="identity")+
    facet_wrap( ~ .data[[facet_factor]],nrow=rows)+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 90), panel.border = element_rect(color = "grey", fill = NA, size = 1))
}

f_pval <- function(x){
  pchisq(-2*sum(log(x)),df = 2*length(x),lower.tail = FALSE)
}
scale_ab <- function(x,a,b){
  (x-min(x))/(max(x)-min(x))*(b-a)+a
}

#define Fisher's method for combining p-values
f_pval <- function(x){
  pchisq(-2*sum(log(x)),df = 2*length(x),lower.tail = FALSE)
}
scale_ab <- function(x,a,b){
  (x-min(x))/(max(x)-min(x))*(b-a)+a
  }
