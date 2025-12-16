function [response] = firstOrderFilter(prevResponse, command, dt, tau)
    eVal = exp(-dt/tau);
    response = eVal*prevResponse + (1-eVal)*command;
end