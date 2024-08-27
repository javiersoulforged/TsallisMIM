
Renyi.empirical <- function(y){

freqs.empirical = function(y)  return( y/sum(y) )

unit="log"

entropy.plugin = function(freqs, unit=c("log", "log2", "log10"))
{
   unit = match.arg(unit)

   freqs = freqs/sum(freqs) # just to make sure ...

   H = -log(sum( ifelse(freqs > 0, freqs^2, 0) ))

   if (unit == "log2")  H = H/log(2)  # change from log to log2 scale
   if (unit == "log10") H = H/log(10) # change from log to log10 scale

   return(H)
}

return( entropy.plugin(freqs.empirical(y), unit=unit) )

}