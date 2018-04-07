CDF   0   
      time             site_id       SGP    facility_id       C1     sds_mode      b1     missing_value         -9999      
sample_int        0.12 seconds   averaging_int         
30 minutes     instruments       0Anemometer:Fill Windmaster Pro, IRGA: LICOR 7500   time_stamp_details        RThe time stamp currently reflects the _start_ of the 30 minute sampling interval.      boom_direction        K0 degrees; wdir is accurate to within an offset of approximately 10 degrees    sign_convention      fc_corr, h, le, mean_g_soil are positive away from surface, r_net is positive towards the surface. That is r_net = h+le+g. Also, r_net = (r_down_short_hemisp+r_down_long_hemisp) - (r_up_short_hemisp+r_up_long_hemisp). fc_corr is defind such that (fc_corr-fc_wpl_h-fc_wpl_le) = fc_uncorr     contact_name      %Marc Fischer, LBNL, mlfischer@lbl.gov      PLEASE_READ_THIS     �For all publications and presentations, please acknowledge: 'U.S. Department of Energy as part of the Atmospheric Radiation Measurement Program.' The automatic inclusion of a data originator as a co-author is not insisted upon in the ARM program, but the source of any data should be clearly recognized either as a co-author or through an appropriate acknowledgment.  The ARM/LBNL Carbon Project contact for this data set is Marc Fischer (mlfischer@lbl.gov). Also please note that we may make adjustments to the data to incorporate adjustments to calibration scales or other issues.  Users should contact Marc Fischer (mlfischer@lbl.gov) to inquire about planned data releases and revisions. Whenever possible, we would appreciate receiving preprints for publications that use the data to insure that the quality and limitations of the data are accurately represented. Your questions and comments are welcome.      MENTOR_QC_FIELD_INFORMATION      �The qc flags use binary, or bitmapped, values (ie values of 0,1,2,4,8,etc.) to note each relevant condition. These can be summed for more detailed qc information. The qc values are:
QC VALUE:  EXPLANATION:
 0-         value not suspect
 1-         value missing
 2-         value below minimum or above maximum or value is +-infinity
 4-         one or more dependencies failed; see the 'dependency' attribute of variable; if a dependency is a 
            spike count, it fails if it is > 100, if a dependency is a variable, it fails when the variable's qc flag
            is not equal to 0
 8-         value has large variance; see 'large variance condition' attribute of variable
16-         value susbect because it has more than 100 spikes; see 'dependency' attribute of variable;
            currently this flag applies to t,q,c
32-         value suspect because of another condition, see 'special condition' attribute of variable.
            Currently only applies to fc_corr, ustar.   	arm_field               irga_serial       C     sonic_serial      B     nrlite_serial         B     history       �created by code 4msonicb1tob1met.c, version r11, with operating system RedHat Linux, kernel 2.4.18-18.7.x, i686 on Apr 09 2015, 08:56:00 GMT      k   	base_time                string        29-Oct-2012,00:00:00 GMT       	long_name         Base time in Epoch     units         #seconds since 1970-1-1 0:00:00 0:00         M�   time_offset                 	long_name         Time offset from base_time     units         'seconds since 29-Oct-2012,00:00:00 GMT          M�   yyyydddhhmmss                   units         yyyydddhhmmss      	long_name         start of integration interval           M�   doy                 units         fractional days    	long_name         fractional day of the year          M�   fc_corr                 units         umol m-2 s-1   	long_name         WPL corrected CO2 flux     	valid_min         ����   	valid_max               
dependency        !nspk_unrot_w,mean_t,mean_q,mean_c           M�   
qc_fc_corr                  units         unitless   	long_name         qc flag for fc_corr         M�   fc_wpl_h                units         umol m-2 s-1   	long_name         #additive WPL H correction to c flux    	valid_min         ����   	valid_max                    M�   qc_fc_wpl_h                 units         unitless   	long_name         qc flag for fc_wpl_h        M�   	fc_wpl_le                   units         umol m-2 s-1   	long_name         $additive WPL LE correction to c flux   	valid_min         ����   	valid_max                    N    qc_fc_wpl_le                units         umol m-2 s-1   	long_name         $additive WPL LE correction to c flux        N   h                   units         W m-2      	long_name         corrected sensible heat flux   	valid_min         ���8   	valid_max               
dependency        nspk_unrot_w,mean_t         N   qc_h                units         unitless   	long_name         qc flag for h           N   le                  units         W m-2      	long_name         WPL corrected latent heat flux     	valid_min         ����   	valid_max               
dependency        nspk_unrot_w,mean_t         N   qc_le                   units         unitless   	long_name         qc flag for le          N   
mean_rot_u                  units         m s-1      	long_name         mean horizontal wind speed          N   qc_mean_rot_u                   units         unitless   	long_name         qc flag for mean_r          N   mean_t                  units         degree C   	long_name         5mean sonic temperature (t), i.e. virtual temperature       	valid_min         ����   	valid_max            2   large_variance_condition      abs(mean_t)/sqrt(var_t) < 2         N    	qc_mean_t                   units         unitless   	long_name         qc flag for mean_t          N$   mean_q                  units         mmol m-3   	long_name         mean H2O density   	valid_min         <#�
   	valid_max         D���   large_variance_condition      abs(mean_q)/sqrt(var_q) < 2         N(   	qc_mean_q                   units         unitless   	long_name         qc flag for mean_q          N,   mean_c                  units         mmol m-3   	long_name         mean CO2 density   	valid_min         A (�   	valid_max         A��        N0   	qc_mean_c                   units         unitless   	long_name         qc flag for mean_c          N4   mean_p                  units         kPa    	long_name         mean IRGA pressure     	valid_min            ^   	valid_max            e        N8   	qc_mean_p                   units         unitless   	long_name         qc flag for mean_p          N<   	var_rot_u                   units         (m s-1)2   	long_name         variance of u      	valid_min                	valid_max            d        N@   qc_var_rot_u                units         unitless   	long_name         qc flag for var_r           ND   	var_rot_v                   units         (m s-1)2   	long_name         variance of v      	valid_min                	valid_max            d        NH   qc_var_rot_v                units         unitless   	long_name         qc flag for var_r           NL   	var_rot_w                   units         (m s-1)2   	long_name         variance of w      	valid_min                	valid_max                    NP   qc_var_rot_w                units         unitless   	long_name         qc flag for var_r           NT   var_t                   units         C2     	long_name         variance of t      	valid_min                	valid_max                    NX   qc_var_t                units         unitless   	long_name         qc flag for var_t           N\   var_q                   units         (mmol m-3)2    	long_name         variance of q      	valid_min                	valid_max           p        N`   qc_var_q                units         unitless   	long_name         qc flag for var_q           Nd   var_c                   units         (umol m-3)2    	long_name         variance of c      	valid_min                	valid_max         >���        Nh   qc_var_c                units         unitless   	long_name         qc flag for var_c           Nl   wdir                units         degrees    	long_name         horizontal wind direction      	valid_min                	valid_max           h        Np   qc_wdir                 units         unitless   	long_name         qc flag for wdir        Nt   theta                   units         degrees    	long_name         rotation to zero w     	valid_min                	valid_max           h        Nx   qc_theta                units         unitless   	long_name         qc flag for theta           N|   phi                 units         degrees    	long_name         rotation to zero v     	valid_min                	valid_max           h        N�   qc_phi                  units         unitless   	long_name         qc flag for phi         N�   ustar                   units         m s-1      	long_name         friction velocity      	valid_min                	valid_max            
        N�   qc_ustar                units         unitless   	long_name         qc flag for ustar           N�   Lmoni                   units         meters     	long_name         Monin-Obukhov length scale     	valid_min         ����   	valid_max           #(        N�   qc_Lmoni                units         unitless   	long_name         qc flag for Lmoni           N�   bar_pres                units         kPa    	long_name         barometric pressure    	valid_min            ^   	valid_max            e        N�   qc_bar_pres                 units         unitless   	long_name         qc flag for bar_pres        N�   t_air_upper                 units         degree C   	long_name         !mean air temperature, upper level      	valid_min         ����   	valid_max            2        N�   qc_t_air_upper                  units         unitless   	long_name         qc flag for t_air_upper         N�   t_air_lower                 units         degree C   	long_name         !mean air temperature, lower level      	valid_min         ����   	valid_max            2        N�   qc_t_air_lower                  units         unitless   	long_name         qc flag for t_air_lower         N�   rh_upper                units         %      	long_name         #mean relative humidity, upper level    	valid_min                	valid_max            d        N�   qc_rh_upper                 units         unitless   	long_name         qc flag for rh_upper        N�   rh_lower                units         %      	long_name         #mean relative humidity, lower level    	valid_min                	valid_max            d        N�   qc_rh_lower                 units         unitless   	long_name         qc flag for rh_lower            N�   z_upper              units         meters     	long_name         height of upper RH/T sensor         M�   z_lower              units         meters     	long_name         height of lower RH/T sensor         M�   mean_g_soil                 units         W m-2      	long_name         -average of up to four soil heat flux sensors       	valid_min         ����   	valid_max            �   special_condition         3less than two out of four sensors read valid values         N�   qc_mean_g_soil                  units         unitless   	long_name         qc flag on mean_g_soil          N�   stderr_mean_g_soil                  units         none   	long_name         standard error of mean_g_soil           N�   z_g_soil             units         meters     	long_name         depth of soil heat flux probes          M�   r_net                   units         W m-2      	long_name         net solar radiation (NR-lite)      	valid_min         ����   	valid_max                    N�   qc_r_net                units         unitless   	long_name         qc flag for r_net           N�   r_tot                   units         W m-2      	long_name         5total downwelling shortwave radiation (LI-200 sensor)      	valid_min         ����   	valid_max           �        N�   qc_r_tot                units         unitless   	long_name         qc flag for r_tot           N�   r_down_short_hemisp                 units         W m-2      	long_name         -downwelling shortwave hemispheric irradiance       	valid_min         ����   	valid_max           �        N�   qc_r_down_short_hemisp                  units         unitless   	long_name         qc flag for r_down_short_hemisp         N�   r_up_short_hemisp                   units         W m-2      	long_name         *upwelling shortwave hemispheric irradiance     	valid_min         ����   	valid_max           �        N�   qc_r_up_short_hemisp                units         unitless   	long_name         qc flag for r_up_short_hemisp           N�   r_down_long_hemisp                  units         W m-2      	long_name         <downwelling longwave hemispheric irradiance; t_rad corrected   	valid_min         ����   	valid_max               
dependency        t_rad           N�   qc_r_down_long_hemisp                   units         unitless   	long_name         qc flag for r_down_long_hemisp          N�   r_up_long_hemisp                units         W m-2      	long_name         :upwelling longwave hemispheric irradiance; t_rad corrected     	valid_min         ����   	valid_max               
dependency        t_rad           N�   qc_r_up_long_hemisp                 units         unitless   	long_name         qc flag for r_up_long_hemisp        N�   t_rad                   units         degree K   	long_name         5radiometer body temperature (for longwave correction)      	valid_min            �   	valid_max           J   
dependency        t_rad           N�   qc_t_rad                units         unitless   	long_name         qc flag for t_rad           O    ppfd                units         umol m-2 s-1   	long_name         (photosynthetic photon flux density (PAR)   	valid_min         ����   	valid_max           �        O   qc_ppfd                 units         unitless   	long_name         qc flag for ppfd        O   mean_t_soil_upper                   units         degree C   	long_name         -average of two upper soil temperature sensors      	valid_min         ����   	valid_max            2   special_condition         (not one of two sensors read valid values        O   qc_mean_t_soil_upper                units         unitless   	long_name         qc flag for mean_t_soil_upper           O   stderr_mean_t_soil_upper                units         none   	long_name         Sstandard error of mean_t_soil_upper; note maximum of two values used in calculation         O   mean_t_soil_middle                  units         degree C   	long_name         .average of two middle soil temperature sensors     	valid_min         ����   	valid_max            2   special_condition         (not one of two sensors read valid values        O   qc_mean_t_soil_middle                   units         unitless   	long_name         qc flag for mean_t_soil_middle          O   stderr_mean_t_soil_middle                   units         none   	long_name         Tstandard error of mean_t_soil_middle; note maximum of two values used in calculation        O    mean_t_soil_lower                   units         degree C   	long_name         -average of two lower soil temperature sensors      	valid_min         ����   	valid_max            2   special_condition         (not one of two sensors read valid values        O$   qc_mean_t_soil_lower                units         unitless   	long_name         qc flag for mean_t_soil_lower           O(   stderr_mean_t_soil_lower                units         none   	long_name         Sstandard error of mean_t_soil_lower; note maximum of two values used in calculation         O,   mean_m_soil_upper                   units         cm3 cm-3   	long_name         +average of four upper soil moisture sensors    	valid_min                	valid_max         ?      special_condition         /less than two of four sensors read valid values         O0   qc_mean_m_soil_upper                units         unitless   	long_name         qc flag for mean_m_soil_upper           O4   stderr_mean_m_soil_upper                units         none   	long_name         $standard error of mean_m_soil_upper         O8   mean_m_soil_lower                   units         cm3 cm-3   	long_name         +average of four lower soil moisture sensors    	valid_min                	valid_max         ?      special_condition         /less than two of four sensors read valid values         O<   qc_mean_m_soil_lower                units         unitless   	long_name         qc flag for mean_m_soil_lower           O@   stderr_mean_m_soil_lower                units         none   	long_name         $standard error of mean_m_soil_lower         OD   precip                  units         
mm 30min-1     	long_name         mean precipitation     	valid_min            ^   	valid_max            e        OH   	qc_precip                   units         unitless   	long_name         qc flag for precip          OL   lat              units         degrees    	long_name         latitude of instrument location         M�   lon              units         degrees    	long_name          longitude of instrument location        M�   alt              units         meters     	long_name         #height of tower base from sea level         M�   ppfd_up                 units         umol m-2 s-1   	long_name         2Upwelling photosynthetic photon flux density (PAR)     	valid_min         ����   	valid_max           �        OP   
qc_ppfd_up                  units         unitless   	long_name         qc flag for ppfd_up         OT   zm               units         	meters         	long_name         #height of instrument from tower bas         M�   zrad             units         meters     	long_name         height of radiation instruments         M�   z_t_soil_upper               units         meters     	long_name         depth of upper t soil probes        M�   z_t_soil_middle              units         meters     	long_name         depth of middle t soil probes           M�   z_t_soil_lower               units         meters     	long_name         depth of lower t soil probes        M�   z_m_soil_upper               units         meters     	long_name         depth of upper m soil probes        M�   z_m_soil_lower               units         meters     	long_name         depth of lower m soil probes        M�P�� @|(�@	���L��Bk����jC�� @���@   �L�;����  ���;L��        B}Hi�  C�� ��{�  ����O  �<U�  ����  �?��'  �@"��  �A��Z  �C���  �A{T?  �BĨd  �=��  �=V9�  �<�JQ  �=�*P  �@���  �:���  �B���  �?{pB  �BT��  �=��  �A�L  �B��+  �@;d  �A*�\  �B�\  �B�  ��'�,    ?^�.�E\)  ��<  ���E�  �>��  �C���  �C��  �C�    ��C��  �A34  �=��A��  �=���A\)  �<u� <UϪ  ��< �<  !��<       ���x�  �@�      B}Hi�׀ C����k�  ���J  ��+b1  ��S�V  ��l�  �@[�  �A�}�  �C�6�  �A|��  �B���  �=��O  �=H�.  �<���  �=��  �>�,2  �:3�  �B�tl  �?2�l  �B]*  �=��Q  �A)4�  �B���  �@I&�  �A�\  �B/�  �B\)  ��=q    ?s���.�
  ��<  ��X�`  �?��  �C�ab  �C�K�  �C�ff  ����T  �A
�\  �=�\ A
    �=� A
��  �<u� <?	�  ��< �<  !��<       ���%  �@�      B}Hi��  C��U>���  ��c�  ����  ���m5  ����  �@ w8  �A���  �C�?�  �A}9^  �B��m  �>#  �>H��  �="�  �=�	w  �?���  �:q  �C v,  �?�h  �BN'R  �> }�  �A�6E  �B�և  �@7�w  �A��  �B/��  �B�R  ����    ?D2��G{  ��<  ��<(�  �@%  �C�+�  �C�>�  �C�33  ���/  �A
Q�  �=���A	�  �>   A
ff  �    <<�  ��< �<  !��<       ���n�  �@�     B}Hi�H� C�� >+�  �����  �:�y�  ���t  ��%�  �@14b  �A�4V  �C���  �A~P  �B��w  �>V�  �>'��  �=896  �=�`  �>���  �9��N  �CW2  �?��  �BF�6  �>=؄  �A�?a  �B�Ձ  �@dZ  �A
�R  �B=�R  �B�\  ���\    ?����[��  ��<  ����  �?�n�  �C��<  �C���  �C��3  ��1'  �@�(�  �=�Q�@��H  �>�@�z�  �    <7,R  ��< �<  !��<       ���-  �@�      B}Hi��  C���>��  ��9i�  ��"m�  ����l  ��k��  �@/�  �A���  �C��h  �A~f  �B���  �>/֢  �=�3  �<�.�  �=���  �?��  �:�q  �Cf�  �@z�  �B2e�  �>U  �Ak�(  �B��N  �?�  �A��  �BKp�  �B��  ���[    ?����V�\  ��<  �����  �@33  �C��l  �C���  �C���  ��ȴ  �@�    �=� @�z�  �>0��@�Q�  �<u� <0J�  ��< �<  !��<       ��(�  �@��     B}Hi��� C��U�0�K  ��O��  �:	N<  ���k;  ��&�  �@�  �A��  �C�zB  �A��  �B��o  �=��  �=���  �<���  �>��  �?;o�  �:�dw  �C�  �@]�  �B#�  �=��B  �Ap  �B��  �?���  �@�
=  �BX
=  �B+��  ����    ?�b��Rz�  ��<  ��|��  �?��  �C�>P  �C���  �C��  ��?}  �@�J  �>��@��
  �>Rn�@��+  �=�<)?)  ��< �<  !��<       ���j  �@�     B}Hi�o  C�� �A�  ��9/  �;�+�  ��~m&  �=��  �@?�  �A4  �C���  �A��]  �B���  �=�W�  �=�H]  �<��   �=��,  �?���  �;~�  �CI�  �@ �R  �B"ي  �=�԰  �Ad?  �B��X  �?��  �@���  �Bc
=  �B3�  �� �:    ?El��>�  ��<  ��Dz�  �?�C�  �C���  �C�H	  �C���  ���  �@���  �>��@��
  �>9� @�    �=�<"	�  ��< �<  !��<       ��p�  �@Ȝ     B}Hi�*� C���<u;�  ��9��  ���p�  ��'�  ��H�  �@O  �A}O�  �C��W  �A��v  �B���  �=�ǉ  �=�+  �<�bi  �=͔Q  �?�  �:ai"  �CqN  �?�;�  �B>:�  �>k�  �AH}�  �B��  �>��  �@ɉ7  �Bh    �B8\)  ��v�    ?v�Z�H  ��<  ��p��  �@6�+  �C��]  �C���  �C��   �����  �@��!  �>�`@��/  �>6��@��H  �<�+ <��  ��< �<  !��<       ��hs  �@�      B}Hi��  C��U�Z��  ����  �<=��  ���lZ  �?�"�  �@lq  �AzV4  �C��  �A��  �B���  �=�l�  �=�X,  �<���  �>
7/  �?�A�  �:��*  �C��  �?ʨ*  �B<t�  �=���  �A���  �B��  �>�  �@��`  �Bo��  �B:��  ����    ?u�2�\ff  ��<  ��F{  �@R��  �C���  �C���  �C�Y�  ��Ƨ�  �@�Z  �=� @�9X  �>C��@���  �<L� <�J  ��< �<  !��<       ���  �@Ϥ     B}Hi��� C�� >��S  ����  �<]a[  ���W  �>���  �@1�  �As��  �C�  �A��  �B���  �>�4  �=�	k  �=:%�  �=��  �@�q\  �:�  �CB4  �@�  �B2�-  �>7�_  �A�m  �B��  ���;d  �@�Q�  �Bt�R  �B8�
  ����    ?���`��  ��<  ��w��  �@(r�  �C�ؔ  �C�A�  �C�33  ���33  �@�(�  �=��@@�~�  �>#T @���  �<T� <��  ��< �<  !��<       ���7L  �@є     B}Hi�Q  C���>϶�  ���9�  ���N  ���Dc  ���c  �@2�S  �Amo  �C�v  �A�av  �B��l  �>��  �=��8  �=8�9  �>c�  �>�y�  �:S�  �C�  �@
�  �B��  �>5<�  �A|p�  �B��  �����  �@�V  �B~��  �B>{  ����    ?�%��_Q�  ��<  �����  �?���  �C�G�  �C�f�  �C�    �����  �@y�  �>��@v5@  �>/�@y�^  �=,@<��  ��< �<  !��<       ���    �@�V     B}Hi�� C��U�� q  ��G�u  ����6  �����  ��G�B  �@1)  �AgL:  �C�S�  �A���  �B���  �>J  �>
�4  �=2  �>"ny  �?��  �;$)-  �C�t  �?��A  �B�2  �>�  �Ar
�  �B�J  ���hs  �@�?}  �B��q  �BG�H  ��    ?����[�  ��<  ����  �?�j  �C��  �C��  �C��f  ���I�  �@FJ  �>��@C�<  �>;h@G  �=O�;���  ��< �<  !��<       ��{�m  �@�     B}Hi��  C�� ��vg  ���*1  ���Ǣ  ��*�  ���C�  �@ ��  �A`X�  �C�Pi  �A�y�  �B���  �=�S�  �=��  �<��  �=׋;  �@a>  �:�   �CD�  �?�|�  �A���  �=���  �A��  �B�i  ��Z  �@M�  �B�    �BV�H  ���7L    ?o���S�\  ��<  ��]?}  �?У�  �C��i  �C�Z  �C��  ��.{  �@��  �>8��@dZ  �>9X@�  �=���<  !��< �<  !��<       ����  �@��     B}Hi�}� C������G  ��'��  �:ސ  ��b��  ����  �@ B�  �A^��  �C�۷  �A��  �B��l  �=�v�  �=��,  �<��&  �=rD�  �?�p�  �9н�  �C1`  �?��o  �A�u  �=�:z  �A<q  �B�\  ��+�
  �@=`B  �B�33  �B\=q  ��߅    ?n'��Q�  ��<  ��V  �@�  �C�C�  �C��!  �C��   ���E�  �?���  �>5<?�t�  �>C�?�G�  �: <�  ��< �<  !��<       ��(��  �@؜     B}Hi�3  C��U��i�  ��ӿp  �:F��  ���  ���ʓ  �@��  �AX��  �C��9  �A��V  �B���  �=�ԗ  �=Y�e  �<En9  �=�,;  �?�C�  �:1  �CW�  �?�r�  �A�C  �=���  �@Ӎs  �B�+  ��E��  �@��  �B���  �Bb�
  ���7M    ?f*)�L    ��<  �� ��  �@M�  �C�	�  �C�d;  �C��3  ��P�`  �?�+  �>6ȴ?�M�  �>J=p?���  �=B�P<V�E  ��< �<  !��<       ����  �@�^     B}Hi�� C�� ��  ���y[  �:�A�  �� </  ��zg  �@��  �A[�s  �C��Z  �A��  �B���  �=���  �=G�  �<�n  �=�!\  �?�@  �:'d�  �C
�  �?���  �B�  �=�u3  �A(P  �B��  ��D�  �@&��  �B���  �Ba�R  ��%    ?T���H�  ��<  ���"�  �@n�R  �C�-z  �C�e�  �C��   ��#�
  �?�7L  �>/?�x�  �>Y�?��F  �=49P<i��  ��< �<  !��<       �����  �@�      B}Hiä  C����-2�  ����  ���/�  ���-  ��ba�  �@�<  �AX�t  �C�Q�  �A�QK  �B���  �=��$  �=�`  �<5F  �=T�?  �>�~�  �9��-  �Cl  �?��  �Aߜ�  �=�Qg  �@�i�  �B�  ��E�T  �@�H  �B���  �Bf�  ���|�    ?/���>�R  ��<  ���V  �@X  �C��(  �C���  �C��   ��#�
  �>�M�  �>|j~>�S�  �>�=p>�n�  �=}�<_�2  ��< �<  !��<       �����  �@��     B}Hi�_� C��U��D  ���0  �:�|�  ����d  ����;  �@��  �AT��  �C��^  �A��  �B���  �=D��  �=��A  �;�-�  �=m�l  �?��}  �:��  �C�N  �?�q  �A�ō  �=���  �@�>�  �B��  ��]O�  �?���  �B���  �Bl
=  ���P    ?Ke��?�  ��<  ��	��  �@Ĝ  �C�+  �C�W  �C�s3  ���1  �>�r�  �>/�>Ձ  �>;d[>�
=  �=	7H<Bzd  ��< �<  !��<       ��'�  �@ߤ     B}Hi�  C�� ���Z  ���wJ  �;� �  ��q  �>eT�  �@%��  �AT��  �C�aM  �A��+  �B���  �=u��  �=��^  �<Fo�  �=v��  �@�9e  �:���  �C#^�  �?׾�  �A�	�  �=�H�  �@�<�  �B��  ��dz�  �?�!  �B�ff  �Bm��  ��ׅ    ?X.��F�  ��<  ���r�  �@<�  �C�{:  �C�$�  �C�L�  ��T��  �?2�   �>bN?/|�  �=��?<j  �<t�<+�  ��< �<  !��<       ��'�  �@�     B}Hi�Ѐ C���>�*j  ��-A�  ���6&  ��iy�  ��#�  �@8�y  �AP�  �C�O�  �A��  �B���  �=�'�  �=�o�  �<��'  �=�V  �@�9  �:�2Z  �C��  �?�  �A�s  �=�:,  �A1H�  �B���  ��i��  �?�n�  �B�    �Bq    ���hs    ?ta��K��  ��<  ����  �@�7  �C��B  �C���  �C�Y�  ���o  �?=O�  �>%`B?6ȴ  �>)��?Ct�  �=B�X<)7  ��< �<  !��<       ��49X  �@�     B}HiȆ  C��U�R�u  ���0�  ����  ���1�  ���P*  �@*r�  �ASzC  �C�  �A��0  �B���  �>p��  �>��x  �=��  �=���  �?��`  �:��k  �Ck�  �?�b  �BP�  �>(��  �A�?I  �B��R  ��ax�  �?�|�  �B�ff  �Bm�  ���`    ?z���K�H  ��<  �����  �@R�  �C��  �C��  �C�ff  ��ě�  �?��  �=��?��D  �>	�\?��^  �<�9`<--�  ��< �<  !��<       ���7L  �@�u     B}Hi�A� C�� �*�  ��Ϋ{  �;��S  ��C  �=���  �@l"  �AOQ�  �C��g  �A�S�  �B���  �=�qq  �=�[�  �<��  �=�hl  �?���  �:X�f  �B�M  �?f۱  �Bee�  �=�Q�  �@��  �B�  ��z�  �?�5?  �B�    �Bw  ��"��    ?gd�?\)  ��<  ��EO�  �@r�  �C��Y  �C��  �C�@   ��49X  �>�;d  �=���>��^  �=�x�>��7  �;�t�<+�]  ��< �<  !��<       ���j  �@�V     B}Hi��  C���>w �  ��8  ���g�  ��xP^  ��v�  �@7�  �AS�q  �C�ZI  �A���  �B��-  �=Ě�  �>�Z  �<�[�  �=�5/  �@���  �:,�  �B�Ml  �?�W[  �B]e(  �=⌶  �@���  �B�  ��`Ĝ  �?�9X  �B�    �Btff  ��z�    ?c1�D(�  ��<  ���z�  �@��T  �C�:  �C���  �C�@   �����  �?Gl�  �=��?<��  �>�8?H�:  �;�<*�  ��< �<  !��<       ��ix�  �@�7     B}Hi˲� C��U��9  ���t�  ���F�  ���P#  ��PQ�  �?���  �AS�?  �C���  �A�D�  �B��8  �=LOl  �=B�  �;��f  �=�c  �At"  �:���  �B�r�  �?1��  �B��N  �<۪�  �?]Ng  �B��  ��hb  �?�    �B�33  �B{�
  ��'O�    ?�9��/�  ��<  ��BJ  �?��m  �C�$i  �C�ל  �C�33  �;o  ��O\)  �>;e��ƨ  �>2-�L��  �;�p<&v!  ��< �<  !��<       �=�h  �@�     B}Hi�h  C�� �<�"  ��#�r  �:\�  ��\&[  ��3�U  �?�3�  �AU�  �C���  �A�?�  �B��  �>*
8  �>~�Y  �;�m�  �=��!  �@���  �:e%�  �C�J  �?�);  �B��  �=3U�  �@F�  �B�
�  ��a%  �?��  �B���  �B}p�  ����    >�0`�2{  ��<  ���Q�  �@l(�  �C��k  �C�Ξ  �C�    �<�  ����  �>�X�0Ĝ  �>��ٿ	�^  �=u<&v!  ��< �<  !��<       �=���  �@��     B}Hi�#� C�«�~,�  ���H5  �;�:�  ��ڀ  �>�e  �@��  �AK��  �C��  �A���  �B���  �<m��  �>�A  �:��  �=�h  �A�6  �;"U�  �B��x  �>S  �B�8�  �<�H  �@ܬ�  �B��  ��{�F  �?��  �B�    �B���  ��(��    ?��?�$�H  ��<  �?��`  �@/�;  �C�  �C�!j  �C��  �@��
  ��9X  �=��x�5@  �=�|��  �:Ĝ <*�  ��< �<  !��<       �?��m  �@��     B}Hi��  C��U�)bz  ���}/  �;b�  ����	  �>>  �@-  �AC؆  �C��  �A�&  �B�/M  �=P��  �=��  �;���  �=�  �@캖  �:%!�  �B�Ol  �?!�  �B���  �=ii}  �@4�  �B�-  ��n{  �?��  �B�    �B�k�  ��    A�Z���  ��<  �B1  �A��H  �C���  �C���  �C��  �B�    �?֗�  �=:^0?�p�  �>�P?�G�  �<� <)��  ��< �<  !��<       �AW�  �@�     B}HiД� C�� �<;  �?&)@  �=Yl�  �A^�|  �@�[g  �@h�  �AF   �C��  �A���  �B�1=  �>R��  �=�  �<�֞  �=��=  �Al�_  �;:`  �B︄  �?�={  �Bp��  �>��  ���/�  �B�?}  ��2-  �@G�  �B�    �Bx�  �B��    B��kA�H  ��<  �CL�  �A�(�  �C���  �C��M  �C��  �C1��  �@���  �>8� @��\  �>��@��y  �=��@�<  !��< �<  !��<       �B(�  �@�     B}Hi�J  C�ʫ��   �@%�J  �=��  �Ba&�  �A4%  �@+[5  �AT�K  �C��e  �A�?6  �B�1�  �>yK  �>���  �=RS  �>�&�  �A
��  �;���  �C	j�  �@��  �B*T�  �>5��  ��Qy  �B�U�  ��U  �@���  �B�33  �B_  �B�ff    A�  Bq
=  ��<  �Ch33  �BM\)  �C��I  �C�M  �C��   �Cצf  �@�G�  �>�@��H  �>\(�@�\  �    �<  !��< �<  !��<       �B{�R  �@�}     B}Hi�� C��U@@p[  �@H23  �>%r  �B�w[  �A��`  �@pW�  �An��  �C��v  �A��F  �B�1�  �?$T�  �>�a  �=�z�  �>���  �A��  �;*�L  �C��  �@@)  �AÓG  �>�6�  ��=��  �B�LJ  �?���  �@�=q  �B��  �BB=q  �A�p�    B�p�C �  ��<  �C��   �B���  �C�b  �C��  �C��3  �Ds3  �A/33  �?5A#�  �>�A'�
  �>���<���  ��< �<  !��<       �B�ff  �@�^     B}HiԻ  C�� �<  ��<  ��<  ��<  ��<  ��<  ��<  ��<  ��<  ��<  ��<  ��<  ��<  ��<  ��<  ��<  ��<  ��<  ��<  ��<  ��<  �B�D  �@1  �A��  �B_
=  �B'�\  �Bi��    B���C;��  ��<  �C��f  �B�33  �C�`�  �C�CH  �C�    �DB    �AFz�  �?Q�A=
>  �>޸PAJz�  �=����<  !��< �<  !��<       �Bי�  �@�?     B}Hi�v� C�ҫ@_d�  �@���  �=kO�  �C-o�  �AX�  �@��P  �A�J  �C��  �A|O�  �B��  �?l��  �?��  �>-޸  �?�  �A��o  �;c��  �C1��  �@):�  �@P�  �>�W  ��<  �B�4:  �@��/  �A:{  �B)�H  �B    ��< B  �< C�Y�  ��<  �C���  �B�    �C�y�  �Cƴ�  �C�&f  �Dq@   �An�\  �?Q�Ak�  �?330Am  �>��;  ��< �<  !��<       �B�ff  �@�      B}Hi�,  C��U? !�  �A	�  �>+��  �CBO  �A؊  �@ˌ�  �A��  �C�0  �Az �  �B���  �?~�J  �?�Ol  �>?��  �?$,�  �A�0�  �:�  �C7�  �@��  �@�u  �>��  �����  �B�!�  �@�I�  �AVff  �A�  �A�  ��< B  �< C�@   ��<  �D9�  �B�ff  �C�ST  �C�5�  �C�33  �D��   �A�(�  �>� A�Q�  �?8Q�A�G�  �?��\<��  ��< �<  !��<       �C�3  �@�     B}Hi�� C�� ?7s[  �A��  �>
��  �CZ`�  �A�Ru  �@�o�  �A���  �C��  �Ay3�  �B��c  �?���  �?��  �>H�  �?+��  �AJ�(  �;y�  �C@^!  �@`�  �AE�  �>ǃL  ���  �B��  �@ҧ�  �AeG�  �A��  �A��  �C�f    BffC��f  ��<  �D @   �Bԙ�  �C�5�  �Cѵ�  �C���  �D�    �A�{  �>uA��R  �?3@A���  �?/\ <�Ȋ  ��< �<  !��<       �C	�   �@��     B}Hiٝ  C�ګ?]��  �A/�D  �>�U  �Cy;   �A�/  �@�S�  �A��  �C�
�  �Ax��  �B��  �?̕�  �?�|<  �>Z��  �?Cӽ  �A+�  �;��  �CFs  �@#��  �A��  �>ӛc  ��˺�  �B��  �@�
=  �Ar=q  �A�  �A�(�  �B���    A���C�ff  ��<  �D+3  �B�33  �C�_k  �CՈ`  �C�33  �D��   �A�ff  �<�� A�    �>�� A��  �>\@<���  ��< �<  !��<       �C�  �@��     B}Hi�X� C��U>��  �A5&�  �> ��  �C�S�  �A�g;  �@� Z  �A�L  �C���  �Aw��  �B�Ç  �?��  �@�  �>PI�  �?dֱ  �@�%�  �;y�g  �CAp  �@(M^  �AQG  �>��0  ���q�  �B��  �@�p�  �A�(�  �A���  �A��\  �B��f    A���C���  ��<  �D1�   �B噚  �C��  �Cٳ  �C��3  �D��   �A��\  �?N�A�G�  �>��A�\)  �>B��<�ҝ  ��< �<  !��<       �C�   �@�     B}Hi�  C�� @���  �ADV�  �=�*�  �C�F  �A���  �@�.�  �A��:  �C�,:  �AvY  �BĝA  �?溕  �?��  �>rxl  �?v�y  �@��  �;Qz�  �C>�  �@.o  �A +$  �>�U�  ���jq  �B��+  �A�H  �A�33  �A�    �Aׅ  �Aw��    BD�C�&f  ��<  �D5�   �B�ff  �C� X  �C�;9  �C�&f  �D�`   �A���  �?<) A�Q�  �?1�A��\  �>Q�<��  ��< �<  !��<       �C    �@�B�    B}Hi�ɀ C��?�ǋ  �A%�h  �=�x.  �CoU�  �A��  �@�Ǵ  �A��  �C���  �At�<  �BćO  �@v�  �@ �.  �>\��  �?f�  �@�9  �;+�  �C?ֆ  �@��  �A=hR  �>�%  ����  �BĬ�  �A(�  �A��
  �A��H  �A���  �A��
    @���C��f  ��<  �D6@   �B�ff  �C��{  �C�%\  �C��   �D�`   �A��  �?b�`A���  �?<(�A�G�  �>���<��p  ��< �<  !��<       �C�  �@�     B}Hi�  C��U?��  �A/`�  �=�g  �C~r�  �A��y  �@��d  �A�t�  �C�h�  �As��  �B�j�  �@�_  �@X�  �>l@  �?Lj�  �@`W  �;��  �CF�&  �@��  �A�2  �>��*  �����  �BĊ�  �A     �A�=q  �A��
  �A�=q  �@�X    @�ƨC�s3  ��<  �D-�3  �B�    �C�^�  �C�i+  �C��   �D�    �A��\  �?Y��A��R  �?]p�A�(�  �>�� <��  ��< �<  !��<       �C �   �@�#�    B}Hi�:� C�� ?��  �A2&�  �=�&�  �C�}3  �A���  �@�s�  �A���  �C�*7  �As1]  �B�>�  �@  �?�*  �>[:!  �?G�  �@s��  �:�L`  �CE��  �@��  �A���  �>�z1  ���2�  �B�o�  �A$��  �A���  �A�Q�  �A��H  ��5�    AS�C��  ��<  �D#�3  �B֙�  �C�ۑ  �C���  �C���  �D��   �A�p�  �?!G�A�33  �?z�@A��  �>�@<���  ��< �<  !��<       �B�33  �@�     B}Hi��  C��@���  �A��  �=�  �C^r�  �AX�a  �@ȸ�  �A�q�  �C�*
  �Ar�  �B�3  �?�c�  �?�TA  �>V�j  �?*��  �@\��  �:�Pw  �CKeI  �@!  �A�*I  �>Ž�  �����  �B�\�  �A(z�  �A�  �A��  �AĸR  ����    A�`�C��   ��<  �D    �Bə�  �C�%w  �C�<�  �C��  �D��   �A��  �?1�A��  �?�38A�=p  �>��<�D�  ��< �<  !��<       �B�ff  �@��    B}Hiᫀ C��U?W߱  �@���  �=|	  �C4��  �A@��  �@�[�  �A��  �C�P#  �Aro�  �B�   �?�~�  �?�*  �>EGW  �?�}  �@1.q  �:�a  �CEl+  �@$^{  �A�aY  �>��+  ���E  �B�F%  �A+33  �A�  �A�33  �A�  ����    A���Cy��  ��<  �C��   �B���  �C��  �C��  �C��3  �Db@   �A�
>  �>޸`A��H  �>uA�
>  �>�Q�<�Ex  ��< �<  !��<       �B���  �@�u     B}Hi�a  C�� ?�:  �@�\  �<���  �B��U  �@��  �@�(c  �A�"�  �C��  �Art�  �B��  �?�u  �@3  �>:,V  �>Zr�  �@1��  �:6�  �CJ�f  �@}R  �A��-  �>�[�  ��q�  �B�:�  �A*=q  �A�  �A��  �A�p�  ��eG�    AE�B�33  ��<  �C�ff  �BH�H  �C�6I  �CӁ�  �C���  �D��  �A�=p  �>z�A���  �>��A�=q  �>��<ݬ�  ��< �<  !��<       �BY��  �@��    B}Hi�� C��>K��  �?�T>  �;�C�  �B$  �@%�  �@˖�  �A��(  �C��F  �Ar�_  �B�C  �?u��  �?'+�  �>�  �=<�C  �@�x
  �9?N�  �CD/U  �@ީ  �A�z�  �>�ߌ  �»3�  �B�/  �A$��  �A��  �A�ff  �A�    ���z�    @���B+p�  ��<  �C3�   �A��  �C�'0  �C�ͷ  �C�    �C�    �A�=q  �=8R A��  �<�� A�{  �>
=�<��  ��< �<  !��<       �B    �@�V     B}Hi��  C��U>���  �?��x  �<�ȼ  �A��  �@T�!  �@�z�  �A��K  �C���  �Ar�K  �B���  �?��  �?�Q  �>
��  �=.�V  �?�S  �9��  �C@dF  �@�  �AFDY  �>�<(  �� +�  �B�  �A#\)  �A�G�  �A���  �A�G�  ��S��    @���A�\)  ��<  �C    �A�(�  �C���  �C��<  �C��  �C���  �A�ff  �=녀A���  �;�� A���  �>(��<�n  ��< �<  !��<       �A�\  �@�ƀ    B}Hi捀 C�� =f��  �?h  �<�Bc  �AU.�  �@fm�  �@�5n  �A�{  �C�S4  �Ar/	  �B��  �?z	�  �?&�D  �> =�  �<��  �@Du�  �8���  �C<��  �@9�  �A��  �>���  ��C+  �B�	�  �A%  �A�    �A�
=  �A�\)  ����    AZ2�A�=q  ��<  �C�3  �A�{  �C��  �C���  �C��  �Cw��  �A�
>  �>��A�z  �=#� A�(�  �>�\ <�!m  ��< �<  !��<       �A��  �@�7     B}Hi�C  C�����  ��[ua  �:��   �����  ��@�(  �@���  �A�E0  �C�?�  �Ar�8  �B��~  �?S�o  �>w��  �=�"  �=�Nq  �>��  �9��B  �C-%�  �@}�  �@�Jr  �>�E\  �B�N�  �B��  �A�
  �A�(�  �A�{  �A��  ��U    @K�'��
  ��<  �B'\)  �AG�  �C�9L  �C��a  �C���  �B�    �A�\)  �=�� A��  �=8R A�Q�  �<�� <�`l  ��< �<  !��<       �@�  �@���    B}Hi��� C��U�7��  ��1ɒ  ��_6  ����  ���Y  �@X�/  �A�i  �C��N  �As��  �B��"  �>)#�  �=�  �<��n  �=�~�  �?o�  �:*3�  �C7  �?���  �@EA4  �>
  �A[�~  �B�
�  �A�  �A
=  �B�  �A�  ���(�    ?����i�R  ��<  ��7K�  �?`Ĝ  �C�L�  �C��P  �C��f  �@�;d  �Ap�\  �=aH Ap��  �=L� Aq
>  �=8Q�<���  ��< �<  !��<       �>ȴ9  �