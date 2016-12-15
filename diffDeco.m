%%%%%differential decoder%%%%%
function out = diffDeco(out)
    a1=out(2:1:end);

    a2=out(1:1:end-1);

    out=a1.*a2;

    out(1:2:end)=-out(1:2:end);