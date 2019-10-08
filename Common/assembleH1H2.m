function [H1, H2] = assembleH1H2(df_dx, df_du, dp_dx, dp_du, prevH1,prevH2)
    if (isempty(prevH1) && isempty(prevH2))
        H1 = df_du*dp_dx+dp_du;
        H2 = df_du*df_dx;
    else
        H1 = [prevH1+prevH2*dp_dx;
            df_du*dp_dx+dp_du];
        H2 = [prevH2*df_dx;
            df_du*df_dx];
    end