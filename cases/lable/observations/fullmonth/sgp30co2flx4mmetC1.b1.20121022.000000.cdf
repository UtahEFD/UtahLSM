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
            Currently only applies to fc_corr, ustar.   	arm_field               irga_serial       C     sonic_serial      B     nrlite_serial         B     history       �created by code 4msonicb1tob1met.c, version r11, with operating system RedHat Linux, kernel 2.4.18-18.7.x, i686 on Apr 09 2015, 08:50:26 GMT      k   	base_time                string        22-Oct-2012,00:00:00 GMT       	long_name         Base time in Epoch     units         #seconds since 1970-1-1 0:00:00 0:00         M�   time_offset                 	long_name         Time offset from base_time     units         'seconds since 22-Oct-2012,00:00:00 GMT          M�   yyyydddhhmmss                   units         yyyydddhhmmss      	long_name         start of integration interval           M�   doy                 units         fractional days    	long_name         fractional day of the year          M�   fc_corr                 units         umol m-2 s-1   	long_name         WPL corrected CO2 flux     	valid_min         ����   	valid_max               
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
qc_ppfd_up                  units         unitless   	long_name         qc flag for ppfd_up         OT   zm               units         	meters         	long_name         #height of instrument from tower bas         M�   zrad             units         meters     	long_name         height of radiation instruments         M�   z_t_soil_upper               units         meters     	long_name         depth of upper t soil probes        M�   z_t_soil_middle              units         meters     	long_name         depth of middle t soil probes           M�   z_t_soil_lower               units         meters     	long_name         depth of lower t soil probes        M�   z_m_soil_upper               units         meters     	long_name         depth of upper m soil probes        M�   z_m_soil_lower               units         meters     	long_name         depth of lower m soil probes        M�P���@|(�@	���L��Bk����jC�� @���@   �L�;����  ���;L��        B}Hc   C�  >��W  ��§  ��WI  ����  ����  �@��9  �B��  �D�?)  �A^/  �B�J�  �?�2#  �?&�9  �>fg�  �=c�V  �B&��  �9l"�  �C$\�  �@Ce  �Az5�  �>���  �CF  �B�Q�  �A��  �A���  �B=�H  �Ba33  �� r�    ?����SG�  ��<  ��i�^  �<��  �C�ĕ  �C�y  �C��3  �=�%  �A�  �=�\�A߸R  �=� A���  �;�� =A�7  ��< �<  !��<       ���l�  �@�      B}Hcۀ C��?TH  ���jn  �=$�  ��/`  ���  �ABN  �B�  �D�]�  �A]�V  �B�Q�  �@kG�  �?�$�  �>��6  �=�  �A��  �96xm  �C.Ш  �@��  �@���  �?�:  �C���  �B�W
  �A�z�  �A��  �B?��  �Bc33  ��:��    ?����d�\  ��<  ��!hs  �?yX  �C��  �C��  �C���  ����  �A߅  �=#� Aߏ\  �=8R A�z�  �<u� =[��  ��< �<  !��<       ���`B  �@�      B}Hc�  C�U>�:J  ���fD  �<��  ���   ����  �AR�  �B�W  �D��  �A^p�  �B�h�  �@��  �?���  �>�n�  �=W�<  �@֖�  �9���  �C,��  �@��  �@�+  �?Y  �C��d  �B�e�  �A���  �A߮  �BA�  �Bdp�  ��=�    ?�]`�ZQ�  ��<  ��a�^  �      �C��  �C�u  �C���  ��ě�  �A��  �=#� A��  �=#� A�
>  �<#� =h�  ��< �<  !��<       ����
  �@�     B}HcL� C� ��ױ  �����  ��e�  ���S*  ��j�  �A&O  �B ��  �D�Q�  �A_{N  �B�lW  �@y2  �?�W�  �>ɉ�  �=e�  �A���  �:3��  �C)�p  �@#ة  �A&�  �?�  �DG�  �B�t:  �A��\  �A��  �BMp�  �Bq�  ��;��    ?��_�Q
=  ��<  ��=�  �      �C���  �C�]�  �C�s3  ���C�  �A׏\  �=8R Aי�  �=#� A�p�  �<#� =oI�  ��< �<  !��<       ��T��  �@�      B}Hc  C�
�>���  ���x  ���V�  ���-  ����  �Ak�  �A�6�  �D��7  �A`P�  �B�tp  �?�R�  �?X��  �>��3  �=��  �A��z  �91ʗ  �C*�  �@d	  �A�  �?\�  �C�DK  �B��>  �A�=q  �A�p�  �BSz�  �BxQ�  ��5��    ?���KQ�  ��<  ��0�`  �      �C��  �C�}  �C�&f  ���o  �A��  �=#� A��  �=#� A�    �;�� =},|  ��< �<  !��<       ����
  �@��     B}Hc�� C�U?
	�  ����o  �<���  ��Mj  ���   �A��  �A���  �D�e�  �A`e�  �B��V  �@"  �?~)�  �>���  �<��  �@
t�  �9)�  �C,�  �@�  �@��  �?(  �C��7  �B���  �A��H  �A�\)  �BTp�  �By=q  ��4�    ?����O��  ��<  ���  �>2-  �C��/  �C�  �C��f  ���o  �A�{  �=\ A��  �=#� A�    �<u� =�0�  ��< �<  !��<       ��D��  �@�     B}Hc
s  C� �J>	  ����r  �<�ս  ����  ��ʚ  �AO  �A�ps  �D��z  �A`��  �B���  �@d  �?��  �>׽�  �<�]�  �@��  �9i
  �C/.&  �@��  �@�;o  �?\�  �Dԃ  �B���  �A��  �A�ff  �BV��  �B{p�  ��;��    ?��c�Sff  ��<  ��#33  �=��P  �C�;�  �C���  �C�ٚ  ��o  �A�=p  �<�� A�G�  �=\ A��  �<u� =���  ��< �<  !��<       �;D��  �@Ȝ     B}Hc.� C��>��O  ����  ����:  ���  ���FS  �A��  �A��  �D���  �A`��  �B���  �@��  �?��  �>˚>  �<���  �@�)M  �8���  �C0�  �@:  �@F*  �?��  �D�  �B���  �A���  �A�(�  �B\Q�  �B���  ��5��    ?����R\)  ��<  ��*n�  �<e`B  �C��L  �C���  �C�ٚ  ��o  �A�    �<�� A�
>  �=] A��H  �<u� =�X�  ��< �<  !��<       ����
  �@�      B}Hc�  C�U>��  ����E  �<�7�  ��dA  ����  �Avp  �A���  �D�oS  �Aa�  �B���  �?�e  �?j��  �>��  �<�,  �?��  �9
��  �C-�  �@W�  �@�.�  �?�  �C��d  �B��r  �A��  �A��  �B`�R  �B�G�  ��*$�    ?�=x�P
=  ��<  ��$�  �=C�  �C��'  �C��  �C��3  �      �Aˏ\  �=\ A˙�  �=#� A�z�  �<u� =�S  ��< �<  !��<       �      �@Ϥ     B}Hc�� C� >� �  ���'   �=	�  ��]/  ��\Q  �A�  �A���  �D�;�  �Aa�  �B��  �?�c/  �?L��  �>���  �=v�  �@W�'  �9yh  �C*3P  �@�:  �A��  �?�  �C���  �B���  �A�Q�  �A�Q�  �Bb��  �B�u�  ��(�    ?�NH�OG�  ��<  ��     �=  �C��q  �C�_  �C�s3  �      �A��  �=#� A��  �=#� A��
  �<#� =� �  ��< �<  !��<       �      �@є     B}HcU  C����O�  ���D�  ��Cв  ��/�y  ����\  �@��;  �A��O  �D�+�  �Aa%�  �B���  �?�n�  �?GpB  �>�o  �=E��  �@�Nm  �9r%�  �C(   �@p�  �A>�  �?�1  �C�X�  �B��'  �A�{  �A˅  �BbQ�  �B�.  ��)+    ?�W��PG�  ��<  ��%/  �>���  �C��  �Cߪ�  �C��   �      �A���  �=L� A��
  �=8R A�  �<u� =���  ��< �<  !��<       �      �@�V     B}Hc� C�U>��  ���E�  �= �x  ��0��  ���y  �@�כ  �A���  �D�7�  �AaqI  �B���  �?�ܣ  �?W�e  �>�a�  �=+�W  �@�}�  �9-1�  �C+$  �@ܧ  �A��  �?�]  �C��P  �B��F  �A���  �A�    �Bc(�  �B���  ��#�    ?��"�S�  ��<  ��$Z  �?K�  �C�	�  �Cߪ  �C��   �      �A�G�  �=#� A�G�  �=#� A�33  �<#� =�p�  ��< �<  !��<       �      �@�     B}Hc�  C�  ���t  ���?9  �<a�w  ��+��  ���_  �AL�  �A�]!  �D�;  �Aa��  �B��   �?�|�  �?�<  �>�_  �=[;:  �@}�X  �9Ͷ�  �C.��  �@��  �@�$,  �>���  �C�D�  �B��F  �A�    �A�{  �Ba=q  �B��=  ��(Ĝ    ?��_��  ��<  ��@A�  �>��h  �C��  �C߃�  �C��   �      �A�{  �=\ A��  �=#� A�    �<u� =�&�  ��< �<  !��<       �      �@��     B}Hc�� C�"�>Ȅ3  ����  �<Չ.  ��8��  ��S��  �@���  �A��P  �D���  �Ab��  �B���  �?ӕ�  �?>�X  �>���  �=;�  �@Ί�  �9u<  �C.`�  �@�  �@���  �? N�  �C|`&  �B��L  �A��  �Ař�  �Bg�H  �B�    ��&��    ?����f    ��<  ��p��  �=ix�  �C�Ô  �C�s�  �C�&f  �;o  �A¸R  �=#� A�  �=\ A£�  �<#� =�B[  ��< �<  !��<       �      �@؜     B}Hc7  C�%U>�~�  ���Lp  �<�fo  ��4��  ��.��  �@�*_  �A�s�  �D�W�  �Ab�{  �B���  �?��C  �?$=  �>�$5  �=7�N  �AW�=  �9D<�  �C1��  �@	4�  �@�m  �>��  �C]��  �B���  �A��  �A��  �BhQ�  �B�=q  ��(�R    ?�x��p=q  ��<  ��Q%  �>���  �C��y  �C���  �C��  �;o  �A���  �<�� A��H  �=#� A��R  �<#� =ě�  ��< �<  !��<       �      �@�^     B}Hc� C�( >�
�  ���2�  �=���  ��GM�  �@��  �@�q�  �A�N�  �D�o�  �Ab�  �B���  �?�=�  �?��  �>�f�  �=BhH  �A��  �9@�)  �C4?  �@  �>|3  �>�W  �C2i  �B�Ö  �A�  �A�z�  �Bg{  �B���  ��'��    ?�w��v    ��<  ��PĜ  �?3t�  �C��+  �C�
G  �C��   �      �A�    �=\ A��  �=#� A���  �<#� =ě�  ��< �<  !��<       �      �@�      B}Hc�  C�*�>�6  ��	�  �=6�  ��R�  ���:  �@�H�  �A�Ti  �D���  �Ac �  �B��;  �?�B�  �?3w�  �>���  �=Xy�  �AAAv  �9J��  �C4�;  �@��  �?2B�  �>���  �CD�  �B��  �A�
=  �A��
  �Be    �B�ff  ��,�`    ?����y    ��<  ��e�  �>�X  �C���  �C��   �C��   �      �A�G�  �=\ A�Q�  �=#� A�34  �<u� =�,�  ��< �<  !��<       �      �@��     B}Hcc� C�-U>�=�  ���  �=ow�  ��V?  ���2  �A  �A��7  �D��  �Ab��  �B��p  �?�de  �?@�V  �>��  �=_�4  �@*�h  �9Y�  �C9HS  �@e�  �@�
>  �>�|�  �C/>�  �B���  �A�ff  �A�  �Bb��  �B�=q  ��5��    ?�yF�|��  ��<  ��U��  �?�w  �C�i   �C��  �C��   �      �A�(�  �<�� A�=p  �=#� A��  �<�� =�Z�  ��< �<  !��<       �      �@ߤ     B}Hc  C�0 >��:  ����  �=c;�  ��R=�  ����z  �A��  �A�g  �D�6�  �Ac�  �B���  �@Ll  �?c�7  �>���  �=n˜  �?�e  �9�g�  �C:�  �@yO  �@� s  �>���  �CW�  �B���  �A�  �A��\  �Bd\)  �B���  ��:�t    ?�q,8R  ��<  ���&�  �>T��  �C�l�  �C�8�  �C���  ��o  �A��  �<�� A�34  �=#� A�
>  �<�� =���  ��< �<  !��<       �      �@�     B}HcԀ C�2�>��  ��dL  �=0��  ��T�  ���  �A��  �A�;  �D�4�  �AcM�  �B���  �?׍�  �?<�  �>�'�  �=Q��  �?+}�  �9RU  �C?%�  �@!<�  �A2X�  �??/  �C`a�  �B���  �A���  �A��
  �Bh�R  �B�\)  ��2�    ?�{�\)  ��<  ���I�  �      �C��]  �C���  �C��3  �      �A��\  �=#� A���  �=#� A�z�  �<�� =�<6  ��< �<  !��<       �      �@�     B}Hc�  C�5U>�&i  ��FE  �<�O  ��R��  ����  �Aό  �A�r�  �D��z  �Ac�  �B���  �?㽕  �?:��  �>�e�  �=c1I  �?Apw  �9|i�  �C?�)  �@%'  �A>�  �>���  �CU��  �B�Ĝ  �A���  �A��  �Bn{  �B�=q  ��-�m    ?�C��}(�  ��<  ����  �=��
  �C�/�  �C٠  �C��   �      �A��  �=#� A�    �=#� A��
  �<�� =�H�  ��< �<  !��<       �;D��  �@�u     B}HcE� C�8 >��f  ���D}  �;�,g  ��?@�  ���B�  �AǗ  �A�M>  �D���  �Adf  �B��A  �?�Y  �?Bhp  �>���  �=L�*  �?��  �9o�}  �CB�  �@�j  �Aa܅  �>�^�  �C4H�  �B��{  �A�ff  �A�\)  �Br(�  �B�ff  ��,��    ?�d4�z��  ��<  ��ix�  �=��m  �C�Q  �C�^�  �C�s3  �;ě�  �A�Q�  �=#� A�\)  �=8R A�33  �<�� =�K]  ��< �<  !��<       �      �@�V     B}Hc�  C�:�>�  ���9T  �<Q��  ��?�8  ����  �AF  �A桿  �D��y  �Ad\3  �B��  �?���  �?w��  �>�-�  �==��  �@�`l  �9N�1  �CA�0  �@a  �A]
  �>�ݿ  �CZef  �B���  �A�\)  �A�z�  �Bu{  �B�    ��/    ?��d��  ��<  ��/|�  �>hs  �C�!k  �C�{�  �C�L�  �=�\)  �A��R  �=\ A���  �=\ A��  �<�� =��  ��< �<  !��<       �<#�
  �@�7     B}Hc�� C�=U>�|  ����  �<�l�  ��&�0  ��R�  �A  �A��1  �D��!  �Adb�  �B��P  �@ W  �?VG  �>���  �=�  �@�D  �9BS�  �C@(  �@�
  �A@ru  �?Q�  �C��  �B��  �A�{  �A��  �Bv
=  �B���  ��3
>    ?Cpf�Mff  ��<  ��'�;  �=ix�  �C�K�  �C٭�  �C��   �=q��  �A���  �<�� A���  �<�� A��  �<u� =�>C  ��< �<  !��<       �;��
  �@�     B}Hc l  C�@ >�"w  ����  �;,}�  ��B  ���f  �A݃  �A���  �D��  �Ad�  �B��1  �?�A�  �?Uu  �>��.  �= �  �@��C  �9'1�  �C@g  �@��  �AFq�  �>��  �C�h	  �B��  �A�{  �A�  �B{    �B�33  ��.      ?I?�L    ��<  ��/  �<�j  �C�n8  �C�~y  �C�L�  �?��  �A�p�  �<�� A�z�  �<�� A�ff  �<u� =�%F  ��< �<  !��<       �=�7L  �@��     B}Hc!'� C�B�>�/�  ���֕  �<O
�  ����  ��k^�  �Apz  �A�6�  �D�[�  �Ad�m  �B��k  �?�R+  �?Ee�  �>���  �<���  �@b��  �9��  �CDX0  �@!��  �A���  �>�c�  �C���  �B���  �A�(�  �A��
  �B|p�  �B���  ��.r�    ?���7��  ��<  �@z~�  �@r�  �C��9  �C�:  �C�@   �A4��  �A���  �=\ A���  �=\ A���  �<�� >u%  ��< �<  !��<       �?�O�  �@��     B}Hc"�  C�EU>��W  ��r~�  �<��  ����  �>�Q  �AV{  �A�މ  �D�f�  �Ae�  �B��6  �?�*�  �?N&  �>���  �<�5�  �@�rR  �8�P�  �CDcu  �@lF  �A��  �>��1  �C�+�  �B�%  �A�z�  �A��  �B|z�  �B���  ��5V    >�����p�  ��<  �A�=q  �@�Ĝ  �C�o�  �Cڴ�  �C��   �BN�  �A�\)  �<�� A�Q�  �<�� A�\)  �<u� >)_  ��< �<  !��<       �@���  �@�     B}Hc#�� C�H >!��  �����  �=w�  ��.  �@�b  �A&�D  �A��  �D�Z  �Ae  �B���  �@�v  �?�HZ  �>�  �<5��  �A�  �8�e  �CA<�  �@g  �AS�~  �?   �E���  �B�u  �A���  �A��  �Bz    �B�ff  ��B�    ?fݿ���  ��<  �Bm33  �A-�  �C��  �C��  �C���  �B�    �A�G�  �<u� A�(�  �<�� A�Q�  �<#� >  ��< �<  !��<       �A&�\  �@�     B}Hc%N  C�J�?'TR  �?VZ�  �9��y  �A�+�  �@?}�  �A+��  �A��  �D��  �Ae  �B��  �@�
  �?��I  �>���  �<y>�  �@�vs  �8c�I  �C?��  �@��  �A:k,  �>�L�  ���lx  �B�
  �A��\  �A�  �Bx�\  �B�ff  ��M�    @��4Bp�  ��<  �B���  �A�
=  �C��?  �C޵�  �C�ٚ  �Ct��  �A�p�  �    A�33  �<#� A��\  �;�� >��  ��< �<  !��<       �A��
  �@�}     B}Hc&	� C�MU?5�  �@ s�  �=���  �BH��  �AvnL  �A2Z$  �A��  �D�8]  �Ad�  �B��  �@P�  �?��z  �>�>  �<�z�  �?�J�  �8b?�  �C=�H  �@�  �A$l  �?�v  �ý�  �B�"N  �A�z�  �A��  �Bsz�  �B�ff  ��[�l    @q��BT(�  ��<  �C    �A�z�  �C��  �C�,�  �C��  �C��  �A���  �;�� A�z�  �;�� A��  �    >�?  ��< �<  !��<       �Aə�  �@�^     B}Hc'�  C�P >��  �?ώ  �=>�2  �B!�   �A5��  �A'�+  �A�I�  �D���  �Ae�  �B�0�  �@+�  �?��O  �>�_  �<���  �@�>  �8���  �CB�2  �@$��  �AjC  �?
k�  �òJ�  �B�<j  �A�(�  �A��  �Br�\  �B���  ��F�    ?6��B(�  ��<  �Bי�  �A���  �C��  �C�,W  �C��  �Ce�3  �A��H  �    AĮ  �;�� A�    �;�� >`�  ��< �<  !��<       �A���  �@�?     B}Hc(z� C�R�>T��  �?F��  �<ڌX  �A���  �@�t�  �A
^   �A�A�  �D��?  �Ae,  �B�d�  �?೙  �?r#1  �>p9$  �<@�  �A�   �8AW  �CI_C  �@�  �A��  �>�Y  �Ù~�  �B�aH  �A��H  �A�p�  �Bp\)  �B���  ��m�.    ?��A��  ��<  �B���  �AmG�  �CΟ�  �C�yW  �C�&f  �CO�3  �A�    �<u� A�  �;�� A�(�  �;�� >�  ��< �<  !��<       �A��H  �@�      B}Hc*0  C�UU>��  �?���  �=_x*  �A�xg  �A.�  �A�  �A�-�  �D�X�  �Aed|  �B�k�  �?�$�  �?L��  �>{��  �<�Bm  �@+�G  �8c��  �CH!e  �@&O  �A�&  �>��  ��nZ  �B�o  �A�G�  �A�Q�  �Bl=q  �B�Q�  �����    @8_wB>��  ��<  �B陚  �A�=q  �C�1�  �C���  �C�33  �C�   �A�  �<�� Aƅ  �;�� A�    �;�� >��  ��< �<  !��<       �A�  �@�     B}Hc*� C�X ?���  �?���  �<���  �B&�v  �A eR  �A2�  �A섄  �D���  �Ae,p  �B�m�  �@��  �?�a�  �>z�  �=Upm  �Ah�  �9��G  �CHT�  �@T�  �A��>  �>�-�  ��6r�  �B�v�  �A��\  �A�
=  �Bl=q  �B�{  ����    @�:5B�33  ��<  �C�3  �A�{  �C�$�  �C�o[  �C�&f  �C�Y�  �A�    �=\ A��
  �;�� A�ff  �;�� >(  ��< �<  !��<       �A�  �@��     B}Hc,�  C�Z�>�;  �@V�o  �>�f�  �B�ϵ  �Bs,�  �A�  �A�n�  �D��  �Ad@&  �B�l�  �@g|  �?���  �>�Tq  �>+j  �C���  �9�8�  �CP�  �@	��  �A�Z  �>�^  �����  �B�v�  �A�  �A�
=  �BT�\  �Bu�H  ����    A>��C!ff  ��<  �C��3  �B'�
  �C�[�  �C�  �C��f  �D9�  �A�=q  �>   AָR  �=#� A׸R  �<�� >p;  ��< �<  !��<       �BQ��  �@��     B}Hc-\� C�]U?"R   �@�O  �>��A  �B�U  �B� �  �A5�  �A��  �D�o�  �Ad8u  �B�l>  �@Ѿ  �?ű�  �>���  �>?�l  �Cvjk  �:�]  �CWX  �@��  �Bi_  �?�[  ���n]  �B�q'  �A���  �A���  �BA
=  �B_z�  ����    A1�JC3L�  ��<  �C��  �B2=q  �C�!�  �C�%j  �C��   �D,�  �A�Q�  �>.@A݅  �=�\ Aޙ�  �=aG =��j  ��< �<  !��<       �B^    �@�     B}Hc/  C�` @��,  �@��  �@��  �C�  �C���  �@���  �A�*t  �Di#q  �AdY  �B�l\  �@���  �A[�  �>Cֳ  �?�  �F� �9�d�  �Cy	�  �?��^  �B�A  �<-�U  ��@�9  �B�o�  �A��\  �AǮ  �B1�R  �BL�  ���%    @��CBљ�  ��<  �CN33  �A�(�  �C�KA  �C�v  �C��3  �C��   �A�(�  �>ff�A��  �=�
�A��  �=�
 =���  ��< �<  !��<       �B��  �@�B�    B}Hc/̀ C�b�>f��  �?�F�  �>Ph�  �B1y�  �A���  �@m�  �Aޭ�  �DM%k  �Aj��  �B�n�  �?G��  �?W"  �=���  �=�=�  �CETz  �:Au�  �C�*�  �?�Ό  �BĪ�  �>�7=  ��2�  �B�u?  �An�\  �A�(�  �B)p�  �BA��  ��J1    @�82B��  ��<  �C�3  �A�Q�  �C��A  �C���  �C��3  �C�Y�  �AÙ�  �>W
@A�  �=�
�A�\)  �=�H =k�  ��< �<  !��<       �Aң�  �@�     B}Hc1�  C�eU>؁l  �?��  �=�T�  �A��  �A`^H  �@�!� �A��.  �DvW�  �Af.O  �B�o�  �A_��  �@}P  �>o�  �?��f  �F%�� �=�J  �CR��  �@>�  �A�vJ  �>��  ��A@*  �B�u�  �A�{  �A��H  �B>
=  �BZ�
  ��PM�    ?�1A�\)  ��<  �B���  �AJ�R  �C�G�  �C��  �C�Y�  �CGL�  �A��  �=u� Aʙ�  �=\ A��H  �<#� =�k�  ��< �<  !��<       �A��\  �@�#�    B}Hc2>� C�h �W�  �@��{  �?�(�  �B憁  �Cs)�  �@	:A �A�x�  �DTA&  �Aj  �B�l%  �?�?�  �A�  �=�J_  �?���  �Ek�  �<爯  �C��|  �?�6�  �B���  �>��  ���}  �B�e`  �AeG�  �A�z�  �B2\)  �BK=q  ��m�    @FFA�=q  ��<  �B���  �A ��  �CϴT  �C�P\  �C���  �C2�3  �A�  �=�� A�G�  �=#� A���  �=aG =a��  ��< �<  !��<       �Ap��  �@�     B}Hc3�  C�j�?��  �@)�  �=С4  �Ba  �A�TL  �AȘ  �A��5  �Dg�  �Ab2;  �B�O�  �@  �?��B  �>��'  �=d�  �C�d�  �9D�  �C;L�  �@K  �@��  �?�L  ��a��  �B�Su  �A�    �Aң�  �B"�  �B<33  ����v    @\sgB���  ��<  �C:�  �A�G�  �Cѡ  �C�0  �C��3  �C�L�  �A��  �    AظR  �    A�G�  �    =�-�  ��< �<  !��<       �B33  �@��    B}Hc4�� C�mU>�b  �@6�d  �>
�E  �B�2�  �AØ  �A/�  �A��  �D`{  �Aa�C  �B�1  �?�A  �?��Z  �>�i�  �=�+�  �Cd  �9M��  �C:�n  �@��  �@ԍ�  �>��  �� ��  �B�N�  �A��R  �Aҏ\  �B   �B5(�  ���b    @DKLBә�  ��<  �CT�  �A���  �CΞ�  �C��)  �C�&f  �C��3  �Aڣ�  �<#� A�=q  �    A���  �<#� =��Y  ��< �<  !��<       �B    �@�u     B}Hc6e  C�p ?8��  �@���  �=���  �Bܤ�  �A���  �@���  �A��,  �D[�  �Aa��  �B�2�  �?��  �?�^  �>h�  �>VI�  �B�+~  �9�S  �C4��  �@c�  �?D�_  �>�9y  �aD  �B�@   �A��  �Aՙ�  �BG�  �B+��  ��ɺ    A	�C:�   ��<  �C�@   �B`��  �C�o�  �C켐  �C���  �D*�  �A�Q�  �>
=@A�ff  �<�� A�G�  �>��=�  ��< �<  !��<       �Bz33  �@��    B}Hc7 � C�r�>��  �@���  �=��  �C%  �A���  �@�+  �A�b  �DP��  �Aa�  �B�'�  �?  �?m�Y  �>5B2  �>��|  �C$�  �:�D  �C7p�  �@*@�  �@\&�  �>�G�  ��G�  �B�*  �A��
  �A�(�  �B��  �B�  ����    Ad�C=�3  ��<  �Cæf  �B��  �C�8�  �C��Y  �C�33  �D6�   �A�{  �?�A�(�  �<�� A�\  �>�f`=�SP  ��< �<  !��<       �B�ff  �@�V     B}Hc8�  C�uU>�z#  �@�X  �<��q  �Bi�t  �A�&  �@Ԧ�  �A�#�  �DE��  �A`ܨ  �B��  �?��k  �?��  �>BA�  �=��  �C2��  �9,�   �C5�;  �@��  �?�.  �>���  �¦�)  �B�  �A�{  �A�{  �A�G�  �B
  ��g�    A>�tB�33  ��<  �C�&f  �B?�  �C�U�  �C�Sb  �C�@   �C�    �A�(�  �>ٙ�A��
  �>uA�    �>޸@=x��  ��< �<  !��<       �BQ�  �@�ƀ    B}Hc9�� C�x >���  �? k�  ���K  �AKE�  ���V  �@ԉg  �A��}  �DF�d  �A`��  �B�  �?��  �?M�o  �>;�'  �<�y�  �A��$  �8ep�  �C79	  �@5B  �@NB�  �>��.  ���#!  �B��  �A�G�  �Aݮ  �B��  �B  ��\�<    AU%�A�Q�  ��<  �C�   �A�\)  �C���  �C��  �C��f  �C��   �A��  �>���A�z�  �>�� A�{  �>�
@=V  ��< �<  !��<       �A�{  �@�7     B}Hc;G  C�z��Z��  ���/�  �<r��  ��G[�  �?��  �@�>�  �A�N�  �DI��  �Aa�  �B���  �?��  �?�  �=�  �<��4  �Bg�  �9 u�  �C>�  �@��  �A.a  �>�ù  �C]�  �B� �  �A���  �A׮  �B	�
  �BG�  ���%    A	��2�R  ��<  �B�ff  �Aw�
  �C��`  �C��  �C�s3  �B�ff  �AڸR  �>�z�Aٮ  �>�Q�A�
=  �>�  =��  ��< �<  !��<       �AV�\  �@���    B}Hc<� C�}U>G�u  �����  �<�u�  ����1  �> 8  �@��(  �A�(  �DO7�  �AapN  �B�/7  �?��L  �?��  �=��D  �=s�r  �B^��  �9.o  �C6��  �@  �@*'�  �>���  �Bn��  �B�.�  �A�    �A�    �B/�  �B%��  ��4Z    ?����>    ��<  �A	p�  �@    �C��7  �Cߓ�  �C�ٚ  �A�Q�  �A̅  �=�A�p�  �=�
 A̙�  �=\ <�|  ��< �<  !��<       �?޸R  �