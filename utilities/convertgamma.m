function gammaOut = convertgamma(gammaIn,framerate,direction)
switch direction
    case 'toOurs'
        gammaOut = exp(log(gammaIn).*(framerate/1000));
    case 'toTheirs'
        gammaOut = exp(log(gammaIn).*(1000/framerate));
end

end