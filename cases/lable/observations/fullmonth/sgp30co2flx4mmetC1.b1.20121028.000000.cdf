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
            Currently only applies to fc_corr, ustar.   	arm_field               irga_serial       C     sonic_serial      B     nrlite_serial         B     history       �created by code 4msonicb1tob1met.c, version r11, with operating system RedHat Linux, kernel 2.4.18-18.7.x, i686 on Apr 09 2015, 08:55:19 GMT      k   	base_time                string        28-Oct-2012,00:00:00 GMT       	long_name         Base time in Epoch     units         #seconds since 1970-1-1 0:00:00 0:00         M�   time_offset                 	long_name         Time offset from base_time     units         'seconds since 28-Oct-2012,00:00:00 GMT          M�   yyyydddhhmmss                   units         yyyydddhhmmss      	long_name         start of integration interval           M�   doy                 units         fractional days    	long_name         fractional day of the year          M�   fc_corr                 units         umol m-2 s-1   	long_name         WPL corrected CO2 flux     	valid_min         ����   	valid_max               
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
qc_ppfd_up                  units         unitless   	long_name         qc flag for ppfd_up         OT   zm               units         	meters         	long_name         #height of instrument from tower bas         M�   zrad             units         meters     	long_name         height of radiation instruments         M�   z_t_soil_upper               units         meters     	long_name         depth of upper t soil probes        M�   z_t_soil_middle              units         meters     	long_name         depth of middle t soil probes           M�   z_t_soil_lower               units         meters     	long_name         depth of lower t soil probes        M�   z_m_soil_upper               units         meters     	long_name         depth of upper m soil probes        M�   z_m_soil_lower               units         meters     	long_name         depth of lower m soil probes        M�P�u�@|(�@	���L��Bk����jC�� @���@   �L�;����  ���;L��        B}Hh��  C�  �ס�  �<���  �8��k  �?A  �<��g  �?��  �A�1c  �C�Q�  �A~%e  �B��m  �=MAb  �=t�  �:5xc  �=�+  �@��  �;X  �CX�  �?vq
  �Bp  �=�4  ����  �B���  �?�j  �@�ff  �B>ff  �B�  ���C    ?cNx�wp�  ��<  ����`  �;ě�  �C�6�  �C��+  �C��   ��(�  �@f-  �>�\(@^�\  �>�I�@q��  �>�8:ȊH  ��< �<  !��<       ��;"�  �@�      B}Hh��� C���*�T  �=%��  �;i�  �?e�  �>���  �?��r  �A���  �C�;�  �A~��  �B��/  �<�G  �<�^'  �9�M�  �=�0  �@�v  �9�K�  �C|�X  �?*�>  �B� �  �<�1�  ��Yca  �B��L  �? �  �@ț�  �BR�
  �B  ����
    ?Xa�^(�  ��<  ��;ƨ  �@(A�  �C���  �C��	  �C�    ��,�D  �@R�+  �>��@M  �>ܬ@` �  �=��p<ȟ  ��< �<  !��<       ��A%  �@�      B}Hh�i  C�U���y $�>��  �<I�  �@X�@ �?�]" �?� p �A|;k  �C��  �A�]
  �B��I  �=�h�  �=���  �9�Ml  �>��b  �A�ۇ  �;�'�  �C��  �>�-�  �B�:  �=F3  ���o�  �B�  ��|�  �@�9X  �BY��  �B$�\  ����    @8O�L  ��<  ���9  �@o�;  �C�B�  �C�/g  �C�s3  ����  �@�  �?>v�@�  �?Z��@    �=�z�< sX  ��< �<  !��<       ���  �@�     B}Hh�$� C� �L�F $��s��  ��'�  ����� ����� �?�n� �Ap�{  �C�?�  �A���  �B���  �=��  �>¸   �;��-  �>�%�  �B�v�  �;Ւ�  �C�  �?-_�  �C$6%  �=�z  �? *�  �B�i  �����  �@��D  �Bm�R  �B1z�  ����T    @Y��S��  ��<  ��;��  �@^�y  �C�O  �C�yQ  �C�    ��p��  �?�v  �>�$�?��x  �>��N?��  �=["�<T�  ��< �<  !��<       �����  �@�      B}Hh��  C�
�?ffQ  ���Q  ��.�  ��'~�  ���gX  �?���  �Am�  �C�rO  �A�t�  �B�-  �<]�x  �=�p�  �:L�  �>2^�  �A�v�  �<��&  �C�gR  �?6�  �B�F  �<��  �?6��  �B��  ���E�  �@S��  �B�B�  �B@(�  ����7    ?�	S�M    ��<  ��1X  �@b�  �C�}�  �C��  �C��f  ����`  �?ix�  �>�2?mV  �>�`A?���  �<�� �<  !��< �<  !��<       ���;d  �@��     B}Hh��� C�U     $�      �      �     �     �?5� �Ar}�  �C��#  �A��  �B��  �> �  �=		  �      �>��s  �Ah.s  �<��  �C�ڠ  �      �B�j�  �      �F<  �B�  ���  �@_�P  �B~=q  �B=Q�  ����    >��|�C    ��<  �� Ĝ  �@�G�  �C�,  �C��  �C�s3  ���  �����  �>~vȾ��H  �>�����
  �=�V:~G�  ��< �<  !��<       ��T��  �@�     B}Hh�K  C� ��  ��$;N  ��H|�  ��^��  ��A#  �?b� �AmWT  �C��R  �A�9�  �B�.�  �<�r$  �>&g  �;S�  �>O`f  �@���  �;�Jo  �C�o2  �?f��  �C,�c  �<���  �??!�  �B�(�  ��5?  �@N��  �B�p�  �BG    ���-    ?�{��N\)  ��<  ���j  �@�?}  �C��s  �C��  �C�ff  ��T��  �>�  �>��E=�`  �?��>���  �<�;�e�  ��< �<  !��<       ��T��  �@Ȝ     B}Hh�� C��>��� $����  �;.*�  ���� �>pu� �?\X. �Ae�  �C��  �A�f�  �B�1b  �>v�  �=�m  �:E  �?S�  �B��  �<�@�  �A���  �>tğ  �C�  �<��x  �>Z�  �B�0!  ��@ �  �@'|�  �B�ff  �BQ{  ��&    ?�)*�O��  ��<  ��
�  �@�t�  �C���  �C���  �C�&f  ��o  �>J��  �>����|�  �?+dZ>�Ĝ  �=k�;�o  ��< �<  !��<       �����  �@�      B}Hhż  C�U<0x�  ��1�  �;�ʥ  ��oS  �>鰋  �?��  �Ag/  �C��D  �A��p  �B�1  �=B��  �=�-<  �:�;�  �>~�  �@��  �;��Z  �B2��  �?1w�  �CY�  �;��"  �:��  �B�1�  ��"�  �@Bn�  �B���  �BOQ�  ����    ?����K(�  ��<  ��v�+  �@���  �C��v  �C���  �C�&f  �      �>Qht  �>䛦���T  �?#�F>�J  �='�;�=�  ��< �<  !��<       ��T��  �@Ϥ     B}Hh�w� C� ���  ��y�  ���IR  ���l%  ��~�  �?���  �A_�r  �C�*	  �A���  �B�1a  �<�}]  �<�gA  �;$ۿ  �<�.�  �@��h  �:f�b  �B~�   �>R��  �B��  �=�  �@,��  �B�4:  ��E��  �@,z�  �B�33  �BU�
  ��$bN    ?��J{  ��<  ��~5?  �@���  �C�8�  �C��  �C�@   �      �?Q&�  �>�-?
��  �?'
=?g+  �<�;�ܜ  ��< �<  !��<       ���o  �@є     B}Hh�-  C����K0  �����  �9��  ���V  ��;  �?��w  �AZ<W  �C�o�  �A�,c  �B�1�  �<�+  �<��  �;*�u  �<�f�  �@-  �:Bf�  �Be�X  �=Xu  �B�)�  �=4ӯ  �@}�  �B�>v  ��X��  �@��  �B���  �B[(�  ��37L    @�R�K(�  ��<  ��	%  �@���  �C��  �C���  �C�33  ��o  �??}  �>�$>��:  �>�dZ?	�^  �<���;�i  ��< �<  !��<       ���t�  �@�V     B}Hh�� C�U���  ���/  �:-��  ��F�D  �=�t  �?��\  �AW�  �C��  �A�Xf  �B�2�  �<RX  �<V;�  �9��  �<��  �@L  �:�
�  �B`��  �>��  �B��$  �;�W�  �=h�>  �B�A  ��l�  �?�bN  �B���  �Be�  ��)    ?��&�C
=  ��<  ��~�  �@z��  �C�W  �C��u  �C�    ��o  ���$�  �>���ff  �>����G�  �=F��;��  ��< �<  !��<       ��#�
  �@�     B}Hhʞ  C�  �ƽ  ��oI�  ���E9  ���B8  �����  �?��"  �AU�G  �C��r  �A��  �B�1�  �=k��  �<�O�  �:�"  �=��  �?�e�  �:vF�  �B�o�  �>վ  �Bא  �=�  �@L�|  �B�A  ��+  �?�  �B���  �Bk    ��!;d    ?����<ff  ��<  ��9�  �@i��  �C�bI  �C�ȯ  �C��   ��ě�  ���bN  �>�I����w  �?
>��\)  �=�1 ;�q�  ��< �<  !��<       ���j  �@��     B}Hh�Y� C�"����  �<�  �:�  �>E�  �=��]  �?�C  �AU��  �C��f  �A�+  �B�1�  �=>�R  �=�#t  �:&  �=Y\�  �?�;�  �:���  �B�O�  �>��  �B��O  �<�RH  ��<B  �B�<j  ��w�P  �?�O�  �B���  �Bf�  ��+�T    ?,l��8
=  ��<  �����  �@���  �C�{�  �C��  �C���  �      ���V  �??����A�  �?L푿��`  �=�O�;��  ��< �<  !��<       ���1  �@؜     B}Hh�  C�%U��� $�>.�Z  ���  �@hw� ��Lw� �?nY� �AN-/  �C��  �A��x  �B�1�  �=��0  �>���  �:�:�  �>�  �A���  �<���  �B��  �=��  �B��  �=k5s  ���D�  �B�<�  ���  �?bJ  �B���  �Bt��  ���h    ?�ZE�<G�  ��<  ��V  �@F�+  �C  �C�nw  �C�s3  �      ��ٙ�  �?E����w  �?6ff���y  �<���;��8  ��< �<  !��<       �=}�  �@�^     B}Hh�ʀ C�( ��h<  ���M  ��K_7  ���  ��܍�  �?�K  �ANc  �C�*�  �A�pa  �B�1�  �=\#  �<���  �:���  �=.�  �Ao�  �<
O�  �B�   �=��  �B���  �<��X  �?���  �B�8R  ���S�  �?+  �B�    �B�  ����    ?KR��4�  ��<  ���r�  �@���  �C�  �C�Y#  �C�&f  �<�1  ���  �?dZ�'��  �?)�_�K�  �    ;�ߏ  ��< �<  !��<       �>x��  �@�      B}Hhπ  C�*�>$�4  ��=��  �<&��  ��|�|  �?�,  �?݊\  �AF�  �C��  �A��u  �B�1f  �=�[  �=C��  �;��  �=�s�  �AbC�  �;�k  �B���  �>�S  �B�%)  �=Lt$  �@9��  �B�:^  ����  �>�E�  �B���  �B�{  ��>�.    ?��*�=��  ��<  ����  �@�r�  �C  �C�-�  �C��  �<�C�  ���h  �>��^� �  �>����j  �<�C�;���  ��< �<  !��<       �>&�y  �@��     B}Hh�;� C�-U� ��  ���ݫ  �9�=�  �����  ��^��  �?���  �AF�k  �C��x  �A���  �B�1�  �=X��  �<�b�  �;L�!  �<���  �?�LM  �;N��  �B�  �>���  �B���  �<ϻC  �?[�  �B�B  ���~�  �>Õ�  �B���  �B���  ��9x�    ?���5p�  ��<  ����  �@���  �C}�k  �C�J  �C��  �=t�  ����0  �>�dZ��j  �>�����n  �;#� ;�/  ��< �<  !��<       �>�^5  �@ߤ     B}Hh��  C�0 >�Z�  ����;  �:�*�  �����  �����  �?��  �A?�  �C��l  �A��  �B�1p  �<�~q  �=O�/  �<��  �=f*�  �@j�  �<
��  �B�Ř  �=��  �B�:g  �=jL�  �@IO�  �B�>v  �����  ��x��  �B�    �B���  ��6��    @8�@\)  ��<  ��M�-  �@7;d  �C~L  �C�@  �C��  �>Kƨ  �� Q�  �>����  �>]/��K�  �<� �;�@x  ��< �<  !��<       �?7K�  �@�     B}HhҬ� C�2��N7  ��Y��  ��9��  ���mM  �� ��  �?�Ra  �AA�c  �C�^  �A�>t  �B�1�  �==�  �=5��  �:���  �<�r<  �?�;>  �;�[�  �B���  �=�pO  �B�T'  �<��G  �@�@  �B�H�  ���M�  �����  �B�    �B�33  ��
�    ?RC2�5�  ��<  �����  �@b�\  �C|b%  �C�61  �C���  �>��-  ��7+  �>�\(�Ap�  �>����5��  �<49�;�º  ��< �<  !��<       �?k  �@�     B}Hh�b  C�5U��1�  ��N�7  ��q�  ���p�  ���\  �>��, �A:j  �C��  �A��  �B�f�  �>GT7  �>u},  �:�Q/  �>�Uk  �Bc�  �;�X  �B��P  �?��  �B��  �=Ŷ8  �@�$v  �B�Rn  �����  ����`  �B�33  �B�ff  ���ff    ?y��4  ��<  ��;d  �@�w  �C|tq  �C��  �C�s3  �>��  �����  �>�=p����  �>���X  �>���;���  ��< �<  !��<       �?co  �@�u     B}Hh�� C�8 ���  �<�O�  �;(!  �>�l.  �>���  �>j� �A72	  �C�Ȭ  �A���  �B�k�  �<�Hj  �=�  �8�T7  �=S�  �@w��  �:DT�  �C�3$  �?8i�  �C0fI  �<�w  �����  �B�V  ���K�  ����  �BÙ�  �B�    ���7M    ?*(�2Q�  ��<  ���t�  �@]��  �C{U�  �C��  �C��f  �?P��  ����<  �>�/���  �?����^  �=�@;���  ��< �<  !��<       �?���  �@�V     B}Hh��  C�:����� $��њF  ���h$  ��
(j ��L�� �?d� �A.^�  �C�u�  �A���  �B�l�  �=�*�  �=*��  �9���  �>k�  �A{��  �9���  �C�5o  �>�n  �B�ջ  �= ��  �@��  �B�\�  ����H  ���;  �B�ff  �B���  ���{    ?F�x�-�H  ��<  �����  �@��  �Cy�2  �C��  �C�L�  �?n{  ���33  �>�`@��M�  �?V�� �  �=m�@;���  ��< �<  !��<       �?���  �@�7     B}Hh׎� C�=U@ m]  ����  ��yjR  ���՞  ����  �?U'� �A*ٔ  �C�H{  �A��  �B�p�  �=�Ye  �=Nܛ  �:.�  �>�  �AG  �=Q"  �C��X  �?���  �B�Ub  �=m܃  �@4A�  �B�e`  ��陚  ��.V  �B���  �B�    ���ě    ?O���(�R  ��<  ���`B  �@�K�  �Cy��  �C�־  �C�ff  �?WK�  ����j  �>KC����  �>\(���1&  �<�;�y  ��< �<  !��<       �?��j  �@�     B}Hh�D  C�@ ��f  �;��  �:��  �=ډ�  �=�o�  �>��� �A3�_  �C�8d  �A�Q�  �B�o�  �<�J  �=Fv  �8��~  �=W�  �?ҍ�  �<Y�K  �C�i  �>{bv  �C�  �<7v�  ����  �B�a�  ��ڏ\  ���5?  �B���  �B���  �����    ?5��2�  ��<  �      �@��F  �C|��  �C��  �C��f  �?���  ����P  �? �����  �?<I���|�  �=�`@;��A  ��< �<  !��<       �?��  �@��     B}Hh��� C�B�>.�� $���.  �;�5�  ��+�d �>��? �?�\� �A!N
  �C��R  �A�X�  �B�mn  �=U&!  �=\&�  �:3�  �>�π  �AĢg  �=
�g  �B>��  �=�-  �C[H  �<��>  �@�FG  �B�]/  ���ff  ��1x�  �B�    �B�    ����    ?J�\�&  ��<  �@IX  �@��^  �Cz��  �C�Ao  �C��   �A
=  ���b  �>n���Q�  �>vȰ���  �<��;� �  ��< �<  !��<       �@�
  �@��     B}Hh۵  C�EU?��| $���z  ����  ����' ��D� �?��� �A�  �C��  �A���  �B�rn  �?��  �<���  �;]F�  �>D�  �A���  �<�3  �Buh�  �?e��  �B�K�  �=�a  �A~4   �B�j  ���z�  ��3��  �B�33 �B�ff  ��K
>    A�z���
=  ��<  �B8G�  �A�p�  �CzK�  �C��  �C��f  �B�33  ���\(  �=�����  �>C����{  �=���<  !��< �<  !��<       �AeG�  �@�     B}Hh�p� C�H ���  �?�C  �=��  �A�,  �A_��  �=��Q �A��  �C�<�  �A�;^  �BŃ�  �=�8  �>)  �<+Ѿ  �>~�|  �B~�a  �<�ez  �C6�  �A��Z  �@3�  �=�fQ  ��bx�  �B�vF  ���E�  �='�  �B�33  �B��  �B      B�'iATz�  ��<  �CL�  �A���  �C}w  �C�iD  �C���  �C@L�  �@6~�  �?@�j  �?"-@01  �?=O��<  !��< �<  !��<       �B    �@�     B}Hh�&  C�J�AY�  ��V�_  ��Y�  ���m$  ���  �?�� �A$�:  �C���  �A�N  �BŃA  �>ؼ�  �>̔�  �=۴  �?F�.  �C�^  �>"�  �C�f  �?���  �A��  �>-�e  �B���  �B�z^  ��X��  �@5?  �B�    �B��R  �B���    A,��B�  ��<  �Cj��  �BQ�H  �C��  �C�)�  �C��   �C�Y�  �@��D  �>�X@�ff  �?3l@�{  �>z��<  !��< �<  !��<       �Bz�
  �@�}     B}Hh�� C�MU>w�  �@��  �>sO  �B��  �A�Q�  �@�  �A@�.  �C��  �A���  �B�y1  �>�Y�  �>���  �=��  �>�B�  �Ag�1  �<�	  �C �H  �@,�!  �BL��  �>h�j  ����  �Bŀ�  ����
  �@�bN  �B���  �Bz�  ��[G�B   �< B���  ��<  �C��  �B�=q  �C�\�  �C�<�  �C�&f  �DS3  �@�  �?<j�@�/  �>~v�@��  �>��<  !��< �<  !��<       �B�ff  �@�^     B}Hh��  C�P @ �F  �@��7  �=�G/  �B��  �A{b�  �@15  �AT�A  �C�i�  �A�:  �B�q*  �>��  �?�K  �=�˟  �>�qY  �@�nA  �;s  �C	�\  �?�]$  �B(b�  �>y�@  ��G��  �Bł�  ���n�  �@��  �B�    �Bc33  �Bř�B   �< CA    ��<  �C��   �B�33  �C�  �C��  �C�ff  �DH�   �A�  �?�A
�R  �>uAQ�  �>�fp<L��  ��< �<  !��<       �B�    �@�?     B}Hh�R� C�R�@]C�  �@؉�  �=���  �C��  �A��o  �@'-�  �A_�  �C�
�  �A��  �B�r4  �>��  �>�R  �=�,d  �>��  �@ռf  �;(v�  �C��  �@,  �B=	  �>B��  ����  �Bł  �?a%  �@�33  �B|�  �BN�H  �C�    BGffC��f  ��<  �C���  �B�ff  �C���  �C�ԕ  �C�ff  �Dp�   �A3�
  �?��PA0ff  �>�  A-
>  �>�G�<SZ�  ��< �<  !��<       �B���  �@�      B}Hh�  C�UU@�I�  �A:  �={��  �C?��  �Ay*C  �@��  �Aj��  �C��c  �A�SN  �B�t0  �?B^  �>��f  �=�j�  �?G*  �A+�F  �;�(  �CM�  �@FR)  �A��K  �>�*Q  ��;q  �Bŀ   �?�w  �A	��  �Bf��  �B;�H  �B���B   �< C�s3  ��<  �D�   �B�ff  �C���  �C�eu  �C�L�  �D�    �A<�H  �>p��AGp�  �?XQ�AC  �?���<Vwq  ��< �<  !��<       �C33  �@�     B}Hh�À C�X ?�)i  �A!?�  �=z�  �C_�~  �A��G  �@j� �AqZ  �C���  �A~��  �B�s~  �?��a  �?Z�  �>ӎ  �?$&k  �A���  �;��g  �C,��  �@K��  �@���  �=�4�  ����  �BŔ{  �@/�w  �A�H  �BU��  �B+��  �C6ffB   �< C��   ��<  �D l�  �Bϙ�  �C�u�  �C��  �C��  �D�`   �AK�
  �>L��A\z�  �?���Al��  �>�z��<  !��< �<  !��<       �C��  �@��     B}Hh�y  C�Z�?�8�  �A��  �=�3j  �C=  �A�پ  �@��  �Ax��  �C���  �A}�)  �B�n�  �?v�  �?rY:  �>i[  �?(�?  �AD]  �;H�  �C%f  �?�L+  �Ao�  �=�J�  ��i8  �BŃ  �@V�R  �A&�H  �BA�  �B�  �C� B   �< C��3  ��<  �D+�  �B�    �C��,  �C��K  �C��3  �D�`   �An    �?/\0Ai�
  �?/\(Ao
>  �>�p��<  !��< �<  !��<       �C��  �@��     B}Hh�4� C�]U>F�  �A!2  �=��  �Ca��  �A��  �@4�  �A���  �C��A  �A|�k  �B�G]  �?fL+  �?��z  �>�=  �?T�)  �AJy  �;<L'  �C	�L  �@��  �B)��  �>E>l  ��6T  �B�b�  �@��D  �A6=q  �B)z�  �B    �B�33B   �< C��  ��<  �D2@   �B�ff  �C��  �C�<�  �C��  �D�`   �A~�H  �?��A|z�  �>� At��  �>���<  !��< �<  !��<       �CL�  �@�     B}Hh��  C�` ��T  �A=�  �>&[y  �C���  �A��  �@�'  �A��c  �C��^  �A{Ti  �B�1  �?T�  �?��l  �>$��  �?f�U  �A0�  �;7  �C
\  �@V{�  �B7֓  �>j'�  ���,   �B�LJ  �@�j  �AD��  �B�  �Bp�  �B�ff    A�33C��   ��<  �D4�   �Bޙ�  �C�u�  �Cև�  �C���  �D�@   �A�p�  �?��A�z�  �?ҏ\A{34  �?#��<  !��< �<  !��<       �C�  �@�B�    B}Hh襀 C�b�����  �A  �=��  �CY��  �A�6  �@� �A� m  �C�'�  �Az[n  �B�G  �?��h  �?�h�  �>%�:  �?$J  �@Ȝ�  �;�'  �C|  �?ڔ�  �B
�  �=�w  �����  �B�5?  �@�j  �AM�  �BQ�  �B =q  �A݅    @��C�ff  ��<  �D3@   �B�ff  �C���  �C׹�  �C�ٚ  �D��   �A��  �?��A�p�  �?�=pA�\)  �?5�<  !��< �<  !��<       �C    �@�     B}Hh�[  C�eU��V�  �A�{  �=���  �CL�=  �A�<X  �@'1 �A�=�  �C�J�  �Ay�J  �B��O  �?�m  �@C  �>]�  �? aK  �@�  �;��  �C�  �@  �A�3e  �>���  ��2��  �B��  �@�J  �AX    �Bz�  �A�{  �@U    A7�C�@   ��<  �D+Y�  �B�ff  �C��  �C�v�  �C��  �D�    �A�p�  �?�\(A�  �?�\(A��  �?E�<�  ��< �<  !��<       �B���  �@�#�    B}Hh�� C�h ;��(  �@�8U  �=Q^�  �C0�  �AOĻ  �@:�� �A��r  �C�a�  �Ay�  �B��U  �?�G�  �?r�  �>T;  �?(�/  �@�25  �:Œ�  �C�  �@0F  �A��S  �>|\C  ���Ǝ  �B��'  �@��y  �A_33  �Bz�  �A���  ���
=    A�C���  ��<  �D9�  �B�33  �C��'  �C�s<  �C�&f  �D��   �A�33  �?��A�ff  �?�fhA�
>  �?<(�<�SQ  ��< �<  !��<       �B�ff  �@�     B}Hh��  C�j���l�  �@��  �=V�  �C$
  �AKZ3  �@.�S �A�)�  �C�E^  �Ax{K  �Bļ�  �?���  �?��R  �>,�  �?Ť  �@ ~�  �:���  �C��  �@
g�  �A�'  �>X�>  �����  �B��)  �@�j  �Ae��  �B\)  �A�p�  ����    A�g�C��   ��<  �D�   �B�33  �C�0�  �C֛  �C��  �D��   �A���  �?�=pA�Q�  �?��XA�
>  �>8Q�<�  ��< �<  !��<       �B�33  �@��    B}Hh퇀 C�mU���  �@���  �=o�I  �C5FH  �Aa��  �@@�
  �A�d�  �C�l  �Axf�  �Bħ�  �?�E�  �@,*  �>�  �?�;  �?���  �:���  �C�  �@ ��  �B89`  �>�b�  ��'V8  �B��V  �@؛�  �Ai�  �A��  �A�=q  ���#�    A��<C|�f  ��<  �D �3  �B�33  �C��  �C�#B  �C�33  �Dk�   �A�=p  �?\(�A���  �?��A��  �>�  <Ʌ�  ��< �<  !��<       �B�33  �@�u     B}Hh�=  C�p �PLx  �@��  �=��  �B�>  �A��  �@7�M  �A��  �C�p5  �Axv&  �BĞ
  �?c�`  �?��q  �=�jY  �>��  �?9�Q  �:Y��  �C�  �@��  �BG��  �>^�  �����  �B���  �@�r�  �Ai�  �B�  �A�Q�  ��K�    A_ҍCL�  ��<  �C�s3  �Bl�
  �C���  �C�  �C�ٚ  �D#    �A�p�  �?��A��\  �?uA�    �>��<���  ��< �<  !��<       �B}=q  �@��    B}Hh��� C�r���ך  �@zp{  �<��  �B�~"  �@�<  �@B�.  �A��  �C��  �AyU�  �BĚB  �?Jp�  �?�,  �=Ĭ�  �>��\  �?��  �:���  �B���  �?�iC  �BT8g  �>��1  ���H�  �B��J  �@ܬ  �Ak\)  �B33  �A�
=  ���(�    An �C    ��<  �C�s3  �B��  �C���  �C��  �C��   �D$`   �A��H  �>�  A�33  �?Z�@A�ff  �>Q�<�  ��< �<  !��<       �B�.  �@�V     B}Hh�  C�uU?���  �@/_�  �<���  �Bx�  �@�3  �@@�4  �A�W  �C� �  �Ay�h  �BĪ  �?ú  �?��  �=�g�  �>!©  �>�]�  �:N�  �B��   �?�2�  �BP��  �>h��  �����  �B���  �@�(�  �Ai�  �A�33  �A�z�  �>H��    AV�B�
=  ��<  �Cz33  �B4�  �C��  �C�2w  �C��f  �C��3  �A��  �=�� A�(�  �?(��A��
  �>��<�cI  ��< �<  !��<       �B>{  �@�ƀ    B}Hh�i� C�x >ѱ�  �>O!;  �;�3�  �@�o  �?L:t  �@8�S  �A��  �C��  �Ay��  �Bķ�  �>���  �>��1  �=E�  �<��$  �>#�  �9�  �C|  �@%�U  �A�L  �>Z%  ��=�  �B��n  �@�7L  �A`Q�  �A�33  �A�  ��D�    A$30�m�  ��<  �B�    �A��\  �C���  �C���  �C��3  �C;�   �Al��  �>�Av    �>ٙ�As33  �<#� <��
  ��< �<  !��<       �A��H  �@�7     B}Hh�  C�z���5  ��X?  �:��Z  ��X��  ���a   �@,0,  �A��  �C���  �Ay�  �BĘ&  �>�s�  �>�  �=-mW  �=@y�  �>�N�  �9��  �C�   �@��  �B1��  �>)�  �A�v�  �B��  �@�1'  �AQ  �B��  �A��\  ��?��    ?����J�R  ��<  �AŅ  �@�K�  �C��  �C�Tp  �C�L�  �BKG�  �AK��  �=��AL��  �<#� ALz�  �<#� <��o  ��< �<  !��<       �@�E�  �@���    B}Hh�ڀ C�}U�о  ��T��  �<%�  ���Ah  �>�G?  �@-Or  �A�q6  �C���  �Az9  �BĖ�  �>1?  �=�+>  �=D�`  �>�H  �@�ұ  �:8�  �C�0  �@�  �B<+=  �>:��  �A��a  �B���  �@�hs  �A?
=  �B
=  �B{  ��&=p    ?ǸP�T(�  ��<  ���G�  �?�\)  �C�7  �C��  �C���  �@^V  �A4    �=�R A3p�  �>   A4=q  �;�� <l�  ��< �<  !��<       �>`A�  �