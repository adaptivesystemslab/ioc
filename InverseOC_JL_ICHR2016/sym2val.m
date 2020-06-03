function mtxOut = sym2val(mtxIn, feature_use_syms, feature_use_val)
J_curr_subbed_q =          subs(mtxIn,                feature_use_syms.q,    feature_use_val.q);
J_curr_subbed_qdq =        subs(J_curr_subbed_q,      feature_use_syms.dq,   feature_use_val.dq);
J_curr_subbed_qdqddq =     subs(J_curr_subbed_qdq,    feature_use_syms.ddq,  feature_use_val.ddq);
mtxOut =                   subs(J_curr_subbed_qdqddq, feature_use_syms.dddq, feature_use_val.dddq);
